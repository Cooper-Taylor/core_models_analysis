// ModelSEED Reaction-Reversibility Heuristics Explorer

const API = {
  data: (path) => fetch(`/data/${path}`).then((r) => r.json()),
  rxn: (rxnId) => fetch(`/api/rxn/${rxnId}`).then((r) => r.json()),
  panelFba: (body) =>
    fetch(`/api/panel_fba`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    }).then((r) => r.json()),
  reactionImpact: (body) =>
    fetch(`/api/reaction_impact`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    }).then((r) => r.json()),
};

const STATE = {
  manifest: null,
  variantPayloads: {}, // tag -> payload (loaded lazily)
  reactionsPanel: null,
  reactionsOther: null, // lazy
  panel: null,
  panelRxnsets: null,   // {model_id: [seed_rxn,...]} for the 100-model panel
  allModelsRxnsets: null, // {model_id: [seed_rxn,...]} for the full 5,683 DB
  allModelsSummary: null, // {generated_at, n_all_models, variants: {tag: {...}}}
  allModelsVariantRows: {}, // tag -> [{model_id, baseline_flux, variant_flux, delta_flux, ...}]
  rxnModelCounts: null, // {rxn: {panel: n, all: m}} — backs the "models" sort column
  selectedVariant: null,
  selectedRxn: null,
  rxnFilter: 'panel',
  rxnVariantFilter: 'any',     // 'any' / 'none' / specific tag
  rxnFluxImpactedOnly: false,  // subfilter for §3c
  scope: 'panel',              // 'panel' | 'all'
  // ----- panel models tab -----
  panelModels: null,           // {model_id: {organism_name, taxonomy, size stats, ...}}
  panelModelVariants: null,    // {model_id: {tag: {delta_flux, ..., n_changed}}}
  panelPipeline: null,         // {model_id: {tag: {base_flux, singles[], cumulative[]}}}
  panelKeyReactions: null,     // {models: {model_id: {base_flux, reactions[]}}, global: [...]}
  panelGrowthControl: null,    // {models: {model_id: {base_flux, n_essential, reactions[], metabolites[]}}, global, metabolites_global}
  panelSyntheticLethal: null,  // {models: {model_id: {pairs[]}}, global: [...]}
  panelFva: null,              // {models: {model_id: {n_blocked, n_forced, reactions[]}}, global: [...]}
  panelRxnDirEffects: null,    // {models:{mid:{base_flux,reactions[]}}, global, options, option_bounds}
  _rde: null,                  // scratch: currently-rendered reaction-direction-effects table state
  methodCmp: null,             // method_comparison.json (agreement/confusion/dist; KEGG_default-anchored model scopes + wide MSDB scopes)
  methodScope: 'models',       // 'models' | 'models_no_transport' | 'all' | 'no_transport' — reaction scope for the method matrices (defaults to the KEGG_default-anchored core-model view)
  analyticsInit: false,        // analytics tab rendered once
  selectedModel: null,
  pmVariantFilter: 'any',
  pmFluxFilter: 'any',
  pmOverrides: {},             // {rxn: mode} for the per-model live override panel
  baselineMap: null,           // {rxn: dir} from baseline.json (live-only; informational)
  pmCompareTag: null,          // variant whose directions sit beside baseline in the per-model table
  pmCompareDiffs: {},          // {rxn: new_dir} for pmCompareTag
  // updated by fetchHealth() in bootstrap; defaults to true so any race
  // (handler firing before fetchHealth resolves) lands in static-safe behavior.
  staticMode: true,
};

async function fetchHealth() {
  try {
    const r = await fetch('/api/health');
    if (!r.ok) throw new Error('status ' + r.status);
    const j = await r.json();
    STATE.staticMode = (j.static_mode === true);
  } catch (e) {
    console.warn('health probe failed, assuming static mode:', e);
    STATE.staticMode = true;
  }
  // index.html ships with <body class="static-mode">; strip it only on
  // confirmed live mode (never add it from JS — keeps the no-flash UX).
  if (!STATE.staticMode) document.body.classList.remove('static-mode');
}

async function loadAllModelsSummary() {
  if (STATE.allModelsSummary) return STATE.allModelsSummary;
  try {
    const j = await API.data('all_models_variants.json');
    STATE.allModelsSummary = j;
    return j;
  } catch (e) {
    console.warn('all_models_variants.json missing — all-models scope disabled:', e);
    return null;
  }
}

async function loadPanelRxnsets() {
  if (STATE.panelRxnsets) return STATE.panelRxnsets;
  STATE.panelRxnsets = await API.data('panel_rxnsets.json');
  return STATE.panelRxnsets;
}

async function loadAllModelsRxnsets() {
  if (STATE.allModelsRxnsets) return STATE.allModelsRxnsets;
  try {
    STATE.allModelsRxnsets = await API.data('all_models_rxnsets.json');
  } catch (e) {
    STATE.allModelsRxnsets = {};
  }
  return STATE.allModelsRxnsets;
}

async function loadAllModelsVariantFba(tag) {
  if (STATE.allModelsVariantRows[tag]) return STATE.allModelsVariantRows[tag];
  try {
    const j = await API.data(`all_models_variant_fba__${tag}.json`);
    STATE.allModelsVariantRows[tag] = j || [];
  } catch (e) {
    STATE.allModelsVariantRows[tag] = [];
  }
  return STATE.allModelsVariantRows[tag];
}

async function loadRxnModelCounts() {
  if (STATE.rxnModelCounts) return STATE.rxnModelCounts;
  try {
    STATE.rxnModelCounts = await API.data('reaction_model_counts.json');
  } catch (e) {
    STATE.rxnModelCounts = {};
  }
  return STATE.rxnModelCounts;
}

// Number of models (in the given scope) that contain reaction `rxn`.
// Returns null when the count is unavailable for that scope (e.g. the
// all-models counts weren't built), so callers can render an em-dash.
function rxnModelCount(rxn, scope) {
  const e = STATE.rxnModelCounts && STATE.rxnModelCounts[rxn];
  if (!e) return null;
  const v = e[scope];
  return (v === null || v === undefined) ? null : v;
}

// Strip the simple inline HTML markup MetaCyc embeds in pathway names
// (CO<sub>2</sub>, <i>N</i>-acetyl…, etc.) so they don't show as literal
// tags after escaping. We render plain text only.
function stripInlineTags(s) {
  return String(s || '').replace(/<\/?(?:sub|sup|i|b|em|strong)>/gi, '');
}

// Render the pathways field: a list of "Source: pwy (desc); pwy2 (desc2); …"
// strings. Group by source, strip inline tags, and lay out each pathway as a
// compact chip so the (often 40+) entries are scannable instead of one blob.
function renderPathways(pathways) {
  const groups = (pathways || []).map((line) => {
    const idx = String(line).indexOf(': ');
    let src = '', body = String(line);
    if (idx > 0) { src = line.slice(0, idx); body = line.slice(idx + 2); }
    const items = body.split(';').map((x) => stripInlineTags(x).trim()).filter(Boolean);
    return { src, items };
  }).filter((g) => g.items.length);
  if (!groups.length) return '';
  return `<div class="pathways">` + groups.map((g) =>
    `<div class="pw-group">` +
    (g.src ? `<span class="pw-source">${escapeHtml(g.src)}</span>` : '') +
    `<div class="pw-items">` +
    g.items.map((it) => `<span class="pw-item">${escapeHtml(it)}</span>`).join('') +
    `</div></div>`
  ).join('') + `</div>`;
}

// "in panel" prevalence line. Uses reaction_model_counts.json so the panel
// count is out of the 100 descriptive panel models (the old `panel_freq`
// field counts over the ~3,461 growers and was mislabeled "of 100 models").
function renderPrevalence(r) {
  const pc = rxnModelCount(r.id, 'panel');   // of 100 descriptive panel models
  const ac = rxnModelCount(r.id, 'all');     // of 5,683 core models
  const allTxt = (ac != null) ? ` &nbsp;·&nbsp; ${ac.toLocaleString()} of 5,683 core models` : '';
  if (pc != null) {
    return pc > 0
      ? `yes — ${pc} of 100 panel models${allTxt}`
      : `no${ac ? ` — present in ${ac.toLocaleString()} of 5,683 core models` : ''}`;
  }
  // Fallback if reaction_model_counts isn't loaded for this rxn.
  return r.in_panel ? 'yes' : 'no';
}

// -------------------- generic sortable table --------------------
// Render a sortable table into `mount`. Re-renders in place on header click.
//   cols: [{ key, label, numeric?, defaultDir?, thClass?, tdClass?,
//            render(row)->html, sortVal(row)->comparable }]
//   rows: array of row objects
//   state: { key, dir } — mutated in place so sort survives re-renders
//   opts: { limit?, moreNote(hiddenCount)->html, onRender(mount) }
function renderSortableTable(mount, cols, rows, state, opts = {}) {
  if (!cols.some((c) => c.key === state.key)) {
    state.key = cols[0].key;
    state.dir = cols[0].defaultDir || 'asc';
  }
  const col = cols.find((c) => c.key === state.key);
  const dir = state.dir === 'asc' ? 1 : -1;
  const valOf = (c, r) => (c.sortVal ? c.sortVal(r) : r[c.key]);
  const sorted = [...rows].sort((a, b) => {
    let va = valOf(col, a), vb = valOf(col, b);
    // nulls always sort to the bottom regardless of direction
    const na = va === null || va === undefined, nb = vb === null || vb === undefined;
    if (na && nb) return 0;
    if (na) return 1;
    if (nb) return -1;
    if (typeof va === 'string' && typeof vb === 'string') return va.localeCompare(vb) * dir;
    return ((va > vb) - (va < vb)) * dir;
  });
  const limit = opts.limit || sorted.length;
  const shown = sorted.slice(0, limit);
  const arrow = (c) => (c.key === state.key
    ? `<span class="sort-arrow">${state.dir === 'asc' ? '▲' : '▼'}</span>`
    : '<span class="sort-arrow dim">↕</span>');
  const thead = `<thead><tr>${cols.map((c) =>
    `<th class="sortable ${c.thClass || ''}${c.key === state.key ? ' sorted' : ''}" data-key="${escapeHtml(c.key)}">` +
    `${escapeHtml(c.label)} ${arrow(c)}</th>`
  ).join('')}</tr></thead>`;
  const tbody = `<tbody>${shown.map((r) =>
    `<tr>${cols.map((c) =>
      `<td class="${c.tdClass || ''}">${c.render ? c.render(r) : escapeHtml(r[c.key])}</td>`
    ).join('')}</tr>`
  ).join('')}</tbody>`;
  const hidden = sorted.length - shown.length;
  const more = hidden > 0
    ? (opts.moreNote ? opts.moreNote(hidden)
       : `<p class="hint">… ${hidden.toLocaleString()} more not shown.</p>`)
    : '';
  mount.innerHTML = `<table class="changed-by-table sortable-table">${thead}${tbody}</table>${more}`;
  mount.querySelectorAll('th.sortable').forEach((th) =>
    th.addEventListener('click', () => {
      const k = th.dataset.key;
      if (state.key === k) {
        state.dir = state.dir === 'asc' ? 'desc' : 'asc';
      } else {
        state.key = k;
        state.dir = (cols.find((c) => c.key === k) || {}).defaultDir || 'desc';
      }
      renderSortableTable(mount, cols, rows, state, opts);
    }));
  if (opts.onRender) opts.onRender(mount);
}

// Wire any `.collapsible-header` inside `root` to toggle its section, with
// the animation driven by the [data-collapsed] attribute (see style.css).
function bindCollapsibles(root) {
  root.querySelectorAll('.collapsible-header').forEach((hdr) => {
    const sec = hdr.closest('.collapsible');
    if (!sec || hdr.dataset.bound) return;
    hdr.dataset.bound = '1';
    const toggle = () => {
      const collapsed = sec.getAttribute('data-collapsed') === 'true';
      sec.setAttribute('data-collapsed', collapsed ? 'false' : 'true');
      hdr.setAttribute('aria-expanded', collapsed ? 'true' : 'false');
    };
    hdr.addEventListener('click', toggle);
    hdr.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); toggle(); }
    });
  });
}

// -------------------- tab switching --------------------
document.querySelectorAll('nav button').forEach((btn) =>
  btn.addEventListener('click', () => {
    document.querySelectorAll('nav button').forEach((b) => b.classList.remove('active'));
    document.querySelectorAll('.tab').forEach((t) => t.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(`tab-${btn.dataset.tab}`).classList.add('active');
    if (btn.dataset.tab === 'reaction' && !STATE.reactionsPanel) loadReactions();
    if (btn.dataset.tab === 'sandbox' && !STATE.staticMode) initSandbox();
    if (btn.dataset.tab === 'panel-models' && !STATE.panelModels) initPanelModels();
    if (btn.dataset.tab === 'analytics' && !STATE.analyticsInit) initAnalytics();
  })
);

// -------------------- variant browser --------------------
async function loadManifest() {
  if (STATE.manifest) return STATE.manifest;
  const m = await API.data('manifest.json');
  STATE.manifest = m;
  return m;
}

function escapeHtml(s) {
  return String(s || '').replace(/[&<>"']/g, (c) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  }[c]));
}

function revBadge(rev) {
  const map = { '>': 'forward', '<': 'reverse', '=': 'free', '?': 'unknown', 'off': 'off' };
  const cls = { '>': 'good', '<': 'warn', '=': 'accent-2', '?': 'dim', 'off': 'dim' };
  return `<span class="tag-pill" style="color: var(--${cls[rev] || 'text'})">${escapeHtml(rev)} (${map[rev] || rev})</span>`;
}

// Count (model, rxn) pairs for one variant's diff against a given rxnsets map.
function countRxnInstancesImpacted(diffs, rxnsets) {
  const changed = new Set(diffs.map((d) => d.rxn));
  let n = 0;
  for (const rxns of Object.values(rxnsets || {})) {
    for (const r of rxns) if (changed.has(r)) n += 1;
  }
  return n;
}

// Sum rxn-instances impacted for one variant + scope. Returns
// {count, n_models_with_change} computed live from rxnsets, OR taken
// from the all-models summary if scope === 'all' and the summary is loaded
// (no need to load the 8MB rxnsets file).
function rxnInstancesForVariant(v, scope) {
  if (scope === 'all') {
    const blk = STATE.allModelsSummary?.variants?.[v.tag];
    if (blk) {
      return {
        count: blk.n_rxn_instances_touched,
        n_models_with_change: blk.n_models_containing_changed_rxn,
      };
    }
    return null;
  }
  // Panel scope: compute from panel_rxnsets, lazily loaded.
  const rxnsets = STATE.panelRxnsets;
  if (!rxnsets) return null;
  // We need the variant's diff list. The manifest only has summary fields,
  // so use the loaded variant payload (which carries diffs).
  const p = STATE.variantPayloads[v.tag];
  if (!p || !p.diffs) return { count: null, n_models_with_change: null, lazy: true };
  const count = countRxnInstancesImpacted(p.diffs, rxnsets);
  const ch = new Set(p.diffs.map((d) => d.rxn));
  let n_models = 0;
  for (const rxns of Object.values(rxnsets)) {
    if (rxns.some((r) => ch.has(r))) n_models += 1;
  }
  return { count, n_models_with_change: n_models };
}

function statsForVariant(v, scope) {
  // Returns {flux_change, flip} per scope.
  if (scope === 'all') {
    const blk = STATE.allModelsSummary?.variants?.[v.tag];
    if (blk) {
      return { flux_change: blk.n_models_flux_change, flip: blk.n_models_flip };
    }
    return { flux_change: null, flip: null };
  }
  return { flux_change: v.n_models_flux_change, flip: v.n_models_flip };
}

async function renderVariants() {
  const [m, summary, _panelRxnsets] = await Promise.all([
    loadManifest(),
    loadAllModelsSummary(),
    // Awaited (not fire-and-forget) so that the panel-scope rxn-instance
    // column is non-empty on first paint.
    loadPanelRxnsets().catch(() => null),
  ]);

  const nAll = summary?.n_all_models;
  if (nAll) {
    document.getElementById('scope-all-count').textContent = nAll.toLocaleString();
  } else {
    document.getElementById('scope-all-count').textContent = '?';
    document.querySelector('.scope-toggle button[data-scope="all"]').disabled = true;
    document.getElementById('scope-status').textContent =
      'all-models data not built — run scripts/build_all_models_impact.py';
  }

  document.getElementById('manifest-meta').innerHTML =
    ` &nbsp;|&nbsp; ${m.variants.length} variants &nbsp;|&nbsp; ` +
    `${m.n_msdb_rxns.toLocaleString()} MSDB reactions &nbsp;|&nbsp; ` +
    `${m.n_panel_rxns.toLocaleString()} reactions in panel union or changed by ≥1 variant &nbsp;|&nbsp; ` +
    `built ${m.generated_at}`;

  updateScopeColumnLabels();

  // Pre-warm: for panel scope, fetch every variant payload so the diff-based
  // rxn-instances counter is non-lazy. Cheap — the files already exist.
  if (STATE.scope === 'panel') {
    await Promise.all(m.variants
      .filter((v) => v.tag !== 'baseline')
      .map((v) => loadVariant(v.tag).catch(() => null)));
  }

  const tbody = document.querySelector('#variants-table tbody');
  tbody.innerHTML = '';
  m.variants.forEach((v) => {
    const tr = document.createElement('tr');
    tr.dataset.tag = v.tag;
    const inst = v.tag === 'baseline' ? { count: 0, n_models_with_change: 0 }
                                      : rxnInstancesForVariant(v, STATE.scope);
    const st = statsForVariant(v, STATE.scope);
    const fmt = (x) => (x == null ? '—' : x.toLocaleString());
    const fmtInst = (inst && inst.count != null)
      ? `${fmt(inst.count)} <span class="hint" style="font-size:11px" ` +
        `title="${fmt(inst.n_models_with_change)} models in scope contain at least one reaction this variant changes">` +
        `(${fmt(inst.n_models_with_change)} models)</span>`
      : '—';
    tr.innerHTML = `
      <td><span class="tag">${escapeHtml(v.tag)}</span></td>
      <td>${escapeHtml(v.section)}</td>
      <td>${escapeHtml(v.apt_title || v.title)}</td>
      <td class="num">${v.n_changed_vs_baseline.toLocaleString()}</td>
      <td class="num">${fmtInst}</td>
      <td class="num">${fmt(st.flux_change)}</td>
      <td class="num">${fmt(st.flip)}</td>
    `;
    tr.addEventListener('click', () => selectVariant(v.tag, tr));
    tbody.appendChild(tr);
  });
}

function updateScopeColumnLabels() {
  const tag = STATE.scope === 'all' ? 'all models' : 'panel';
  document.getElementById('col-rxn-instances').textContent =
    `rxn instances impacted (${tag})`;
  document.getElementById('col-flux-change').textContent =
    `models flux change (${tag})`;
  document.getElementById('col-grow-flip').textContent =
    `models grow-flip (${tag})`;
}

async function loadVariant(tag) {
  if (STATE.variantPayloads[tag]) return STATE.variantPayloads[tag];
  if (tag === 'baseline') {
    // synthesize a payload for baseline (no diff vs itself)
    const m = await loadManifest();
    const meta = m.variants.find((v) => v.tag === 'baseline');
    const payload = {
      tag: 'baseline',
      title: meta.title,
      apt_title: meta.apt_title,
      description: meta.description,
      citations: meta.citations,
      section: meta.section,
      diffs: [],
      transitions: {},
      n_changed: 0,
      rev_counts: meta.rev_counts,
      panel_fba: [],
      n_models_flip: 0,
      n_models_flux_change: 0,
    };
    STATE.variantPayloads[tag] = payload;
    return payload;
  }
  const p = await API.data(`variants/${tag}.json`);
  STATE.variantPayloads[tag] = p;
  return p;
}

async function selectVariant(tag, tr) {
  document.querySelectorAll('#variants-table tbody tr').forEach((r) => r.classList.remove('selected'));
  if (tr) tr.classList.add('selected');
  STATE.selectedVariant = tag;
  const p = await loadVariant(tag);
  // Per-variant JSON files predate the manifest's apt_title/description/citations
  // fields. Merge them from the manifest (single source of truth) before render.
  const m = await loadManifest();
  const meta = m.variants.find((v) => v.tag === tag) || {};
  if (meta.apt_title && !p.apt_title) p.apt_title = meta.apt_title;
  if (meta.description && !p.description) p.description = meta.description;
  if (meta.citations && !p.citations) p.citations = meta.citations;
  // Per-reaction model counts back the "models" column in the changed-
  // reactions table (and its default sort) — load once, lazily.
  await loadRxnModelCounts();
  // For all-models scope, the per-variant FBA rows are useful — load lazily.
  let allFbaRows = null;
  if (STATE.scope === 'all' && tag !== 'baseline') {
    allFbaRows = await loadAllModelsVariantFba(tag);
  }
  renderVariantDetail(p, allFbaRows);
}

function renderVariantChangeCard(p) {
  const title = p.apt_title || p.title;
  const desc = p.description ? `<p class="vcc-desc">${escapeHtml(p.description)}</p>` : '';
  const cites = (p.citations && p.citations.length)
    ? `<div class="vcc-refs"><span class="vcc-refs-label">References</span>` +
      p.citations.map((c) => `<span class="vcc-cite">${escapeHtml(c)}</span>`).join('') +
      `</div>`
    : '';
  return `
    <div class="variant-change-card">
      <div class="vcc-label">Change</div>
      <div class="vcc-title">${escapeHtml(title)}</div>
      ${desc}
      ${cites}
    </div>`;
}

function renderPanelFluxChart(panel_fba) {
  const sorted = [...(panel_fba || [])].sort((a, b) => b.delta_flux - a.delta_flux);
  if (!sorted.length) {
    return '<p class="hint">No panel FBA data for this variant.</p>';
  }
  const max_abs = Math.max(...sorted.map((r) => Math.abs(r.delta_flux)));
  if (max_abs < 1e-6) {
    return `<p class="hint">All ${sorted.length} panel models have |Δ flux| &lt; 1e-6 under this variant — the heuristic does not change biology in this panel.</p>`;
  }
  const n_pos = sorted.filter((r) => r.delta_flux > 1e-6).length;
  const n_neg = sorted.filter((r) => r.delta_flux < -1e-6).length;
  const n_zero = sorted.length - n_pos - n_neg;
  const n_flip = sorted.filter((r) => r.baseline_grows !== r.variant_grows).length;
  const rows = sorted.map((r) => {
    const flipped = r.baseline_grows !== r.variant_grows;
    const sign = r.delta_flux > 1e-6 ? 'pos' : (r.delta_flux < -1e-6 ? 'neg' : 'zero');
    const halfPct = (Math.abs(r.delta_flux) / max_abs * 50).toFixed(2);
    const flipIcon = flipped
      ? (r.variant_grows
          ? '<span class="pfc-flip-icon up">↑</span>'
          : '<span class="pfc-flip-icon down">↓</span>')
      : '';
    const tooltip =
      `${r.model_id}\n` +
      `baseline: ${r.baseline_flux.toFixed(4)} (${r.baseline_grows ? 'grows' : 'no growth'})\n` +
      `variant:  ${r.variant_flux.toFixed(4)} (${r.variant_grows ? 'grows' : 'no growth'})\n` +
      `Δ growth flux: ${r.delta_flux.toFixed(4)}${flipped ? '   ← grow-status flipped' : ''}`;
    const bar = sign === 'zero'
      ? ''
      : `<span class="pfc-bar ${sign}" style="width:${halfPct}%"></span>`;
    return `<div class="pfc-row${flipped ? ' flipped' : ''}" title="${escapeHtml(tooltip)}">` +
           `<span class="pfc-id">${escapeHtml(r.model_id)}</span>` +
           `<span class="pfc-flip">${flipIcon}</span>` +
           `<span class="pfc-bar-cell">${bar}</span>` +
           `<span class="pfc-val ${sign}">${r.delta_flux >= 0 ? '+' : ''}${r.delta_flux.toFixed(3)}</span>` +
           `</div>`;
  }).join('');
  return `
    <div class="panel-flux-chart">
      <div class="pfc-legend">
        <strong>${sorted.length}</strong> panel models &nbsp;·&nbsp;
        <span class="legend-item"><span class="swatch pos"></span> ${n_pos} gained flux</span> &nbsp;·&nbsp;
        <span class="legend-item"><span class="swatch neg"></span> ${n_neg} lost flux</span> &nbsp;·&nbsp;
        ${n_zero} unchanged &nbsp;·&nbsp;
        ${n_flip} flipped grow-status (↑ became grower, ↓ stopped growing) &nbsp;·&nbsp;
        bar scale: ±${max_abs.toFixed(2)} growth flux units
      </div>
      <div class="pfc-rows">${rows}</div>
    </div>`;
}

function renderAllModelsStatsCard(p, allFbaRows) {
  const summary = STATE.allModelsSummary?.variants?.[p.tag];
  if (!summary) {
    return `<p class="hint">All-models stats unavailable for this variant.</p>`;
  }
  const fmtFlux = (x) => (x == null ? '—' : Number(x).toFixed(4));
  const linkModel = (mid) =>
    mid ? `<code>${escapeHtml(mid)}</code>` : '—';
  const inc = summary.max_increase;
  const dec = summary.max_decrease;
  return `
    <div class="all-models-stats-grid">
      <div class="ams-card">
        <h4>scope</h4>
        <div class="stat">${summary.n_all_models.toLocaleString()}<span class="stat-sub"> models</span></div>
        <div class="stat-sub">${summary.n_models_containing_changed_rxn.toLocaleString()} contain ≥1 changed rxn</div>
      </div>
      <div class="ams-card">
        <h4>growth status</h4>
        <div class="stat">${summary.n_models_grow.toLocaleString()}<span class="stat-sub"> grow under variant</span></div>
        <div class="stat-sub"><span class="diff-down">${summary.n_grew_default_now_not.toLocaleString()}</span> grew in default → don't now</div>
        <div class="stat-sub"><span class="diff-up">${summary.n_not_default_now_grew.toLocaleString()}</span> no growth in default → now do</div>
      </div>
      <div class="ams-card">
        <h4>flux distribution (variant)</h4>
        <div class="stat">mean ${fmtFlux(summary.mean_flux)}</div>
        <div class="stat-sub">median ${fmtFlux(summary.median_flux)}</div>
        <div class="stat-sub">std    ${fmtFlux(summary.std_flux)}</div>
      </div>
      <div class="ams-card ams-card-wide">
        <h4>largest flux increase</h4>
        <div class="stat">${linkModel(inc?.model_id)}</div>
        <div class="stat-sub"><span class="diff-up">+${fmtFlux(inc?.delta_flux)}</span> &nbsp; baseline ${fmtFlux(inc?.baseline_flux)} → variant ${fmtFlux(inc?.variant_flux)}</div>
      </div>
      <div class="ams-card ams-card-wide">
        <h4>largest flux decrease</h4>
        <div class="stat">${linkModel(dec?.model_id)}</div>
        <div class="stat-sub"><span class="diff-down">${fmtFlux(dec?.delta_flux)}</span> &nbsp; baseline ${fmtFlux(dec?.baseline_flux)} → variant ${fmtFlux(dec?.variant_flux)}</div>
      </div>
    </div>
  `;
}

function renderTopMoversTable(rows, maxRows = 100) {
  if (!rows || !rows.length) {
    return '<p class="hint">No models have changed reactions in this variant.</p>';
  }
  const sorted = [...rows].sort((a, b) => Math.abs(b.delta_flux) - Math.abs(a.delta_flux));
  const shown = sorted.slice(0, maxRows);
  const more = sorted.length - shown.length;
  const tbody = shown.map((r) => {
    const sign = r.delta_flux > 1e-6 ? 'pos' : (r.delta_flux < -1e-6 ? 'neg' : 'zero');
    const flipped = r.baseline_grows !== r.variant_grows;
    const flag = flipped ? (r.variant_grows ? ' ↑' : ' ↓') : '';
    return `<tr class="${flipped ? 'flipped' : ''}">
      <td><code>${escapeHtml(r.model_id)}</code></td>
      <td class="num">${Number(r.baseline_flux).toFixed(4)}</td>
      <td class="num">${Number(r.variant_flux).toFixed(4)}${flag}</td>
      <td class="num pfc-val ${sign}">${r.delta_flux >= 0 ? '+' : ''}${Number(r.delta_flux).toFixed(4)}</td>
      <td>${r.baseline_grows ? '✓' : '✗'} → ${r.variant_grows ? '✓' : '✗'}</td>
    </tr>`;
  }).join('');
  const note = more > 0
    ? `<p class="hint">… ${more.toLocaleString()} more models with at least one changed rxn (showing top ${maxRows} by |Δ flux|).</p>`
    : '';
  return `
    <table class="changed-by-table">
      <thead><tr><th>model_id</th><th>baseline growth flux</th><th>variant growth flux</th><th>Δ growth flux</th><th>grow base → variant</th></tr></thead>
      <tbody>${tbody}</tbody>
    </table>${note}`;
}

function renderVariantDetail(p, allFbaRows = null) {
  const pane = document.getElementById('variant-detail');
  let html = `
    <h3>${escapeHtml(p.tag)} — ${escapeHtml(p.apt_title || p.title)}</h3>
    <p class="hint">${escapeHtml(p.section)}</p>
    ${renderVariantChangeCard(p)}
    <h3>Reversibility-count snapshot</h3>
    <div class="transition-grid">
      ${Object.entries(p.rev_counts).filter(([k,_]) => k.startsWith('new_rev')).map(([k, v]) => {
        const rev = k.replace('new_rev=', '');
        return `<div class="t-cell">${revBadge(rev)}<span>${Number(v).toLocaleString()}</span></div>`;
      }).join('')}
    </div>`;

  if (p.tag !== 'baseline') {
    html += `
      <h3>Transitions vs baseline (${p.n_changed.toLocaleString()} rxns changed)</h3>
      <div class="transition-grid transition-grid-wide">
        ${Object.entries(p.transitions || {}).sort().map(([k, v]) => {
          const [from, to] = k.split('->');
          const dir = (to === '>' || to === '<') && (from === '=' || from === '?') ? 'up'
                      : (from === '>' || from === '<') && (to === '=' || to === '?') ? 'down' : '';
          return `<div class="t-cell ${dir}">` +
                 `<span class="t-trans">${revBadge(from)} → ${revBadge(to)}</span>` +
                 `<span class="t-count">${Number(v).toLocaleString()}</span></div>`;
        }).join('')}
      </div>`;

    if (STATE.scope === 'panel') {
      html += `
        <h3>Panel FBA impact</h3>
        <div class="flux-impact-grid">
          <div class="card"><h4>panel size</h4><div class="stat">${p.panel_fba.length}</div></div>
          <div class="card"><h4>models grow-status flipped</h4><div class="stat">${p.n_models_flip}</div></div>
          <div class="card"><h4>models with growth flux Δ &gt; 1e-6</h4><div class="stat">${p.n_models_flux_change}</div></div>
        </div>`;
    } else {
      const sm = STATE.allModelsSummary?.variants?.[p.tag];
      html += `
        <h3>All-models FBA impact <span class="hint">— ${STATE.allModelsSummary?.n_all_models?.toLocaleString() || '?'} core models</span></h3>
        ${renderAllModelsStatsCard(p, allFbaRows)}`;
    }

    const scopeLabel = STATE.scope === 'all' ? 'all models' : 'panel';
    html += `
      <div class="collapsible" id="changed-rxn-collapsible" data-collapsed="false">
        <div class="collapsible-header" role="button" tabindex="0" aria-expanded="true">
          <span class="collapse-caret">▾</span>
          <h3>Top reactions changed
            <span class="hint">— top 50 by models (${escapeHtml(scopeLabel)}); click a column to sort</span>
          </h3>
        </div>
        <div class="collapsible-body"><div class="collapsible-inner">
          <div id="changed-rxn-table-mount"></div>
        </div></div>
      </div>`;

    if (STATE.scope === 'panel') {
      html += `
        <h3>Panel-wide Δ flux <span class="hint">— all ${p.panel_fba.length} models vs baseline, sorted by Δ; hover any row for numbers</span></h3>
        ${renderPanelFluxChart(p.panel_fba)}`;
    } else if (allFbaRows && allFbaRows.length) {
      html += `
        <h3>Top movers across all models <span class="hint">— ranked by |Δ flux|; only models containing ≥1 changed reaction are listed</span></h3>
        ${renderTopMoversTable(allFbaRows)}`;
    }
  }

  pane.innerHTML = html;

  // Animated collapsible sections (the "Top reactions changed" dropdown).
  bindCollapsibles(pane);

  // Sortable "Top reactions changed" table — default sort by model count
  // (descending), in the active scope. Re-binds rxn cross-links on each
  // re-render via onRender.
  const mount = document.getElementById('changed-rxn-table-mount');
  if (mount && p.tag !== 'baseline') {
    const scope = STATE.scope;
    const rows = (p.diffs || []).map((d) => ({
      rxn: d.rxn, base: d.base, new: d.new,
      models: rxnModelCount(d.rxn, scope),
    }));
    const modelsLabel = `models (${scope === 'all' ? 'all' : 'panel'})`;
    const cols = [
      { key: 'rxn', label: 'rxn', defaultDir: 'asc',
        render: (r) => `<a href="#" class="rxn-link" data-rxn="${escapeHtml(r.rxn)}">${escapeHtml(r.rxn)}</a>`,
        sortVal: (r) => r.rxn },
      { key: 'base', label: 'base', defaultDir: 'asc',
        render: (r) => revBadge(r.base), sortVal: (r) => r.base },
      { key: 'new', label: 'new', defaultDir: 'asc',
        render: (r) => revBadge(r.new), sortVal: (r) => r.new },
      { key: 'models', label: modelsLabel, numeric: true, defaultDir: 'desc',
        thClass: 'num', tdClass: 'num',
        render: (r) => (r.models == null ? '—' : r.models.toLocaleString()),
        sortVal: (r) => r.models },
    ];
    const sortState = { key: 'models', dir: 'desc' };
    renderSortableTable(mount, cols, rows, sortState, {
      limit: 50,
      moreNote: (hidden) =>
        `<p class="hint">… ${hidden.toLocaleString()} more not shown. Use the Reaction Explorer to browse.</p>`,
      onRender: bindRxnLinks,
    });
  }

  // Cross-link any remaining rxn IDs in the variant view → reaction explorer
  bindRxnLinks(pane);
}

// Bind every `.rxn-link` inside `root` to jump to the Reaction Explorer.
function bindRxnLinks(root) {
  root.querySelectorAll('.rxn-link').forEach((a) => {
    if (a.dataset.bound) return;
    a.dataset.bound = '1';
    a.addEventListener('click', (e) => {
      e.preventDefault();
      document.querySelector('nav button[data-tab="reaction"]').click();
      setTimeout(() => selectRxn(a.dataset.rxn), 30);
    });
  });
}

// -------------------- reaction explorer --------------------
async function loadReactions() {
  if (!STATE.reactionsPanel) {
    STATE.reactionsPanel = await API.data('reactions_panel.json');
    document.getElementById('rxn-stats').textContent =
      Object.keys(STATE.reactionsPanel).length.toString();
  }
  // Populate the variant filter dropdown (one option per variant)
  const vsel = document.getElementById('rxn-variant-filter');
  if (vsel && vsel.options.length <= 2) {
    const m = await loadManifest();
    m.variants.forEach((v) => {
      if (v.tag === 'baseline') return;
      const o = document.createElement('option');
      o.value = v.tag;
      o.textContent = `${v.tag} — ${v.title}`;
      vsel.appendChild(o);
    });
  }
  // Pre-load panel rxnsets so the flux-impacted-only subfilter works.
  loadPanelRxnsets().catch(() => {});
  // Pre-load per-reaction model counts so the "in panel" prevalence line
  // uses the true 100-panel / 5,683-DB denominators.
  loadRxnModelCounts().catch(() => {});
  renderReactionList();
}

// Cache: tag -> Set of rxn_ids whose direction-change has a downstream
// effect on growth flux in at least one panel model containing it.
const _FLUX_IMPACTED_CACHE = {};
async function fluxImpactedRxnsFor(tag) {
  if (_FLUX_IMPACTED_CACHE[tag]) return _FLUX_IMPACTED_CACHE[tag];
  const p = await loadVariant(tag);
  const rxnsets = await loadPanelRxnsets();
  // Models in panel with |Δ flux| > 1e-6 for this variant
  const movers = new Set((p.panel_fba || [])
    .filter((r) => Math.abs(r.delta_flux) > 1e-6)
    .map((r) => r.model_id));
  // For each changed reaction, check whether ANY panel-mover model contains it.
  const out = new Set();
  const movedRxnSets = {};
  for (const mid of movers) {
    movedRxnSets[mid] = new Set(rxnsets[mid] || []);
  }
  for (const d of (p.diffs || [])) {
    for (const mid of movers) {
      if (movedRxnSets[mid].has(d.rxn)) { out.add(d.rxn); break; }
    }
  }
  _FLUX_IMPACTED_CACHE[tag] = out;
  return out;
}

async function ensureReactionsOther() {
  if (!STATE.reactionsOther) {
    STATE.reactionsOther = await API.data('reactions_other.json');
  }
  return STATE.reactionsOther;
}

async function renderReactionList() {
  const q = document.getElementById('rxn-search').value.toLowerCase().trim();
  const filt = document.getElementById('rxn-filter').value;
  const varFilt = STATE.rxnVariantFilter;
  const fluxOnly = STATE.rxnFluxImpactedOnly;
  const isSpecificVariant = !(varFilt === 'any' || varFilt === 'none');
  const needsOther = (filt !== 'panel' && filt !== 'changed_panel');
  if (needsOther && !STATE.reactionsOther) {
    document.getElementById('reaction-list').innerHTML =
      '<li><em class="hint">Loading reactions index (~4 MB)…</em></li>';
    document.getElementById('rxn-result-count').textContent = '';
    return;
  }
  let entries;
  if (filt === 'panel' || filt === 'changed_panel') {
    entries = Object.values(STATE.reactionsPanel);
  } else {
    const other = STATE.reactionsOther || {};
    entries = [...Object.values(STATE.reactionsPanel), ...Object.values(other)];
  }
  if (filt === 'changed' || filt === 'changed_panel') {
    entries = entries.filter((r) => (r.changed_by || []).length > 0);
  }
  if (isSpecificVariant) {
    entries = entries.filter((r) =>
      (r.changed_by || []).some((c) => c.variant === varFilt));
    if (fluxOnly) {
      const fluxSet = await fluxImpactedRxnsFor(varFilt);
      entries = entries.filter((r) => fluxSet.has(r.id));
    }
  }
  if (q) {
    entries = entries.filter((r) =>
      (r.id || '').toLowerCase().includes(q) ||
      (r.name || '').toLowerCase().includes(q) ||
      (r.definition || '').toLowerCase().includes(q));
  }
  entries.sort((a, b) => (a.id || '').localeCompare(b.id || ''));
  const ul = document.getElementById('reaction-list');
  const extra = isSpecificVariant
    ? ` (variant=${varFilt}${fluxOnly ? ', flux-impacted only' : ''})`
    : '';
  document.getElementById('rxn-result-count').textContent =
    `${entries.length.toLocaleString()} reactions${extra}`;
  ul.innerHTML = '';
  entries.slice(0, 500).forEach((r) => {
    const li = document.createElement('li');
    if (!r.in_panel) li.classList.add('rxn-not-in-panel');
    li.dataset.rxn = r.id;
    li.innerHTML = `<strong>${escapeHtml(r.id)}</strong>${r.changed_by && r.changed_by.length ? ` <span class="tag">Δ${r.changed_by.length}</span>` : ''}` +
                   `<span class="rxn-name">${escapeHtml(r.name || '(no name)')}</span>`;
    li.addEventListener('click', () => selectRxn(r.id));
    ul.appendChild(li);
  });
  if (entries.length > 500) {
    const li = document.createElement('li');
    li.innerHTML = `<em class="hint">… ${entries.length - 500} more (narrow with search)</em>`;
    ul.appendChild(li);
  }
}

async function selectRxn(rxnId) {
  STATE.selectedRxn = rxnId;
  document.querySelectorAll('#reaction-list li').forEach((li) => li.classList.remove('selected'));
  const li = document.querySelector(`#reaction-list li[data-rxn="${rxnId}"]`);
  if (li) li.classList.add('selected');
  const pane = document.getElementById('reaction-detail');
  pane.innerHTML = '<p class="loading">loading…</p>';
  let r = STATE.reactionsPanel[rxnId];
  if (!r) {
    await ensureReactionsOther();
    r = STATE.reactionsOther[rxnId];
  }
  if (!r && !STATE.staticMode) {
    // last-resort server lookup (rxns not in any index)
    try { r = await API.rxn(rxnId); } catch (e) { /* noop */ }
  }
  if (!r || r.error) {
    pane.innerHTML = `<p class="hint">No data for <code>${escapeHtml(rxnId)}</code>.</p>`;
    return;
  }
  // Ensure per-reaction model counts are available for the prevalence line.
  await loadRxnModelCounts();
  renderReactionDetail(r);
}

function renderReactionDetail(r) {
  const pane = document.getElementById('reaction-detail');
  const stoich = (r.stoichiometry || []).map((s) => {
    const cls = s.coef < 0 ? 'coef-neg' : 'coef-pos';
    return `<li class="${cls}">${s.coef >= 0 ? '+' : ''}${s.coef.toFixed(3).replace(/\.?0+$/, '')}  <strong>${escapeHtml(s.cpd)}</strong>@${s.cpt}  ${escapeHtml(s.name || '')}${s.formula ? ' [' + escapeHtml(s.formula) + ']' : ''}</li>`;
  }).join('');

  const changedRows = (r.changed_by || []).map((c) =>
    `<tr><td>${escapeHtml(c.variant)}</td><td>${revBadge(c.base)}</td><td>${revBadge(c.new)}</td></tr>`
  ).join('');

  // Analytic P(direction) badges for the variants we have it for.
  const pdir = r.p_direction || {};
  let pdirHtml = '';
  if (Object.keys(pdir).length) {
    const rows = Object.entries(pdir).map(([variant, p]) => {
      const pct = (x) => `${(100 * x).toFixed(1)}%`;
      const bar = (color, w) => `<div style="background:${color};width:${(100*w).toFixed(1)}%;height:8px;display:inline-block;vertical-align:middle"></div>`;
      return `<tr>
        <td>${escapeHtml(variant)}</td>
        <td class="num">${pct(p.p_forward)}</td>
        <td class="num">${pct(p.p_reverse)}</td>
        <td class="num">${pct(p.p_reversible)}</td>
        <td style="min-width:200px">${bar('var(--good)', p.p_forward)}${bar('var(--accent-2)', p.p_reversible)}${bar('var(--warn)', p.p_reverse)}</td>
      </tr>`;
    }).join('');
    pdirHtml = `<h3>Analytic P(direction) <span class="hint">— §3.6 of the review</span></h3>
      <p class="hint">From the marginal CC normal on ΔG′° (treating concentration term as fixed).
      Green = P(forward), amber = P(reversible), red = P(reverse).</p>
      <table class="changed-by-table">
        <thead><tr><th>variant</th><th>P(fwd)</th><th>P(rev)</th><th>P(reversible)</th><th></th></tr></thead>
        <tbody>${rows}</tbody>
      </table>`;
  }

  // §3d — simplified live sweep: one chosen "changed direction" derived from
  // the selected variant (only variants that actually alter this reaction
  // populate the select). For each panel model that contains this reaction,
  // show its baseline flux and its variant flux side-by-side. No mode picker.
  const changedBy = r.changed_by || [];
  const sweepHtml = STATE.staticMode
    ? `<h3>Per-variant panel sweep <span class="hint">— live FBA</span></h3>
       <p class="hint static-only-notice">Per-variant panel sweep requires live FBA.
       Restart the server with <code>python3 site/serve.py --live</code>
       (requires cobra + the ModelSEEDDatabase + .kbcache).
       See README "Enabling live FBA".</p>`
    : (changedBy.length === 0
      ? `<h3>Per-variant panel sweep <span class="hint">— live FBA</span></h3>
         <p class="hint">No variant changes this reaction's direction, so there is no "changed direction" to compare against the baseline.</p>`
      : `<h3>Per-variant panel sweep <span class="hint">— live FBA</span></h3>
         <p class="hint">Run FBA on every panel model that contains this reaction, comparing
           the cascade baseline against the selected variant's new direction for this reaction.
           First call to a model loads the JSON; expect ~5-20s.</p>
         <div class="run-row">
           <label>variant:
             <select id="rxn-sweep-variant"></select>
           </label>
           <button id="rxn-sweep-run" class="primary">Run sweep</button>
           <span id="rxn-sweep-status" class="hint"></span>
         </div>
         <div id="rxn-sweep-results"></div>`);

  pane.innerHTML = `
    <h3>${escapeHtml(r.id)} <span class="hint">— ${escapeHtml(r.name || '(no name)')}</span></h3>
    <dl>
      <dt>definition</dt><dd>${escapeHtml(r.definition || '')}</dd>
      ${r.equation ? `<dt>equation</dt><dd><code>${escapeHtml(r.equation)}</code></dd>` : ''}
      <dt>is_transport</dt><dd>${r.is_transport ? 'yes' : 'no'}</dd>
      <dt>ΔG′° (kcal/mol)</dt><dd>${r.deltag != null ? r.deltag.toFixed(2) : '—'} ± ${r.deltagerr != null ? r.deltagerr.toFixed(2) : '—'}</dd>
      ${r.ec_numbers && r.ec_numbers.length ? `<dt>EC numbers</dt><dd>${r.ec_numbers.map(escapeHtml).join(', ')}</dd>` : ''}
      ${r.pathways && r.pathways.length ? `<dt>pathways</dt><dd>${renderPathways(r.pathways)}</dd>` : ''}
      ${r.in_panel !== undefined ? `<dt>in panel</dt><dd>${renderPrevalence(r)}</dd>` : ''}
    </dl>
    ${pdirHtml}

    ${stoich ? `<h3>stoichiometry</h3><ul class="stoich-list">${stoich}</ul>` : ''}

    ${changedRows ? `<h3>Changed by ${r.changed_by.length} variant${r.changed_by.length > 1 ? 's' : ''}</h3>
      <table class="changed-by-table">
        <thead><tr><th>variant</th><th>baseline</th><th>variant</th></tr></thead>
        <tbody>${changedRows}</tbody>
      </table>` : `<p class="hint">No variants change this reaction's direction vs baseline.</p>`}

    ${sweepHtml}
  `;

  if (!STATE.staticMode && changedBy.length > 0) {
    // §3d — only list variants that actually change this reaction's direction.
    // Each option carries the new direction; the run handler uses it as the
    // single sweep mode (baseline + that direction = the two flux numbers).
    const sel = document.getElementById('rxn-sweep-variant');
    sel.innerHTML = '';
    changedBy.forEach((c) => {
      const o = document.createElement('option');
      o.value = c.variant;
      o.dataset.newRev = c.new;
      o.textContent = `${c.variant}: ${c.base} → ${c.new}`;
      sel.appendChild(o);
    });

    const REV_TO_MODE = { '>': 'forward', '<': 'reverse', '=': 'free' };
    document.getElementById('rxn-sweep-run').addEventListener('click', async () => {
      const selEl = document.getElementById('rxn-sweep-variant');
      const variant = selEl.value;
      const newRev = selEl.selectedOptions[0]?.dataset.newRev;
      const mode = REV_TO_MODE[newRev] || 'off';   // ? → 'off' as a safe fallback
      const btn = document.getElementById('rxn-sweep-run');
      const status = document.getElementById('rxn-sweep-status');
      const out = document.getElementById('rxn-sweep-results');
      btn.disabled = true;
      status.textContent = 'running FBA sweep…';
      out.innerHTML = '';
      try {
        const t0 = performance.now();
        const res = await API.reactionImpact({
          rxn_id: r.id,
          variant: 'baseline',     // compare against the unaltered cascade...
          modes: [mode],           // ...with the reaction forced to the variant's new direction
        });
        status.textContent = `done in ${(performance.now() - t0) / 1000 | 0}s ` +
          `(${res.n_models} models, ${(res.elapsed_s).toFixed(1)}s server-side)`;
        renderSweepSimple(res, out, { variant, newRev, mode });
      } catch (exc) {
        status.textContent = 'error: ' + exc.message;
      } finally {
        btn.disabled = false;
      }
    });
  }
}

function renderSweepSimple(res, out, opts) {
  // §3d — minimal three-column view: panel model name, default flux,
  // and flux when this reaction is forced to the variant's new direction.
  const { variant, newRev, mode } = opts;
  const modeRows = res.by_mode[mode] || {};
  const mids = Object.keys(res.baseline).sort();
  let nGrewBefore = 0, nGrewAfter = 0, nFlipped = 0, nFluxChanged = 0;
  // Build structured rows for the sortable table.
  const rows = mids.map((mid) => {
    const b = res.baseline[mid];
    const v = modeRows[mid];
    if (b.grows) nGrewBefore += 1;
    if (v && v.grows) nGrewAfter += 1;
    const flipped = v ? (v.grows !== b.grows) : false;
    if (flipped) nFlipped += 1;
    if (v && Math.abs(v.delta_flux) > 1e-6) nFluxChanged += 1;
    return {
      model_id: mid,
      base_flux: b.growth_flux,
      base_grows: b.grows,
      variant_flux: v ? v.growth_flux : null,
      delta_flux: v ? v.delta_flux : null,
      variant_grows: v ? v.grows : null,
      flipped,
      has_variant: !!v,
    };
  });

  out.innerHTML = `
    <div class="flux-impact-grid">
      <div class="card">
        <h4>variant</h4>
        <div class="stat" style="font-size:14px"><code>${escapeHtml(variant)}</code></div>
        <div class="stat-sub">new direction = ${escapeHtml(newRev)} (mode: ${escapeHtml(mode)})</div>
      </div>
      <div class="card">
        <h4>panel growth status</h4>
        <div class="stat">${nGrewBefore} → ${nGrewAfter}</div>
        <div class="stat-sub">${nFlipped} flipped, ${nFluxChanged} growth flux Δ</div>
      </div>
      <div class="card">
        <h4>panel models</h4>
        <div class="stat">${mids.length}</div>
        <div class="stat-sub">containing this reaction</div>
      </div>
    </div>
    <h4>Per-model growth flux <span class="hint">— click a column to sort (asc/desc)</span></h4>
    <div id="sweep-table-mount"></div>
    <p class="hint">"↑" = model became grower under the variant's new direction for this reaction; "↓" = model stopped growing. "✗" = was not a grower in the baseline.</p>
  `;

  const fluxCls = (d) => (d == null ? '' : (d > 0 ? 'diff-up' : (d < 0 ? 'diff-down' : 'diff-zero')));
  const cols = [
    { key: 'model_id', label: 'model_id', defaultDir: 'asc',
      render: (r) => `<code>${escapeHtml(r.model_id)}</code>`, sortVal: (r) => r.model_id },
    { key: 'base_flux', label: 'default growth flux', numeric: true, defaultDir: 'desc',
      thClass: 'num', tdClass: 'num',
      render: (r) => `${r.base_flux.toFixed(4)}${r.base_grows ? '' : ' ✗'}`,
      sortVal: (r) => r.base_flux },
    { key: 'variant_flux', label: `variant growth flux (${newRev})`, numeric: true, defaultDir: 'desc',
      thClass: 'num',
      tdClass: 'num',
      render: (r) => {
        if (!r.has_variant) return '—';
        const flag = r.flipped ? (r.variant_grows ? ' ↑' : ' ↓') : '';
        return `<span class="${fluxCls(r.delta_flux)}">${r.variant_flux.toFixed(4)}${flag}</span>`;
      },
      sortVal: (r) => r.variant_flux },
    { key: 'delta_flux', label: 'Δ growth flux', numeric: true, defaultDir: 'desc',
      thClass: 'num', tdClass: 'num',
      render: (r) => {
        if (!r.has_variant) return '—';
        return `<span class="${fluxCls(r.delta_flux)}">${r.delta_flux >= 0 ? '+' : ''}${r.delta_flux.toFixed(4)}</span>`;
      },
      sortVal: (r) => r.delta_flux },
  ];
  // Default: largest |Δ flux| first — most-impacted models on top.
  const sortState = { key: 'delta_flux', dir: 'desc' };
  renderSortableTable(document.getElementById('sweep-table-mount'), cols, rows, sortState);
}

// Search wiring
document.getElementById('rxn-search').addEventListener('input', () => renderReactionList());
document.getElementById('rxn-filter').addEventListener('change', (e) => {
  STATE.rxnFilter = e.target.value;
  if (STATE.rxnFilter !== 'panel') ensureReactionsOther().then(renderReactionList);
  else renderReactionList();
});
document.getElementById('rxn-variant-filter').addEventListener('change', (e) => {
  STATE.rxnVariantFilter = e.target.value;
  // Hide the flux-impacted-only subfilter entirely until a specific
  // variant is selected — it has no meaning without one, and showing it
  // greyed-out reads as confusing.
  const sub = document.getElementById('rxn-flux-impacted-only');
  const label = document.querySelector('.rxn-flux-filter-label');
  const subEnabled = !(e.target.value === 'any' || e.target.value === 'none');
  label.hidden = !subEnabled;
  if (!subEnabled) sub.checked = false;
  STATE.rxnFluxImpactedOnly = sub.checked;
  renderReactionList();
});
document.getElementById('rxn-flux-impacted-only').addEventListener('change', (e) => {
  STATE.rxnFluxImpactedOnly = e.target.checked;
  renderReactionList();
});

// ===== Scope toggle (panel vs all-models) wiring (Variant Browser) =====
document.querySelectorAll('.scope-toggle button').forEach((btn) =>
  btn.addEventListener('click', () => {
    if (btn.disabled) return;
    const scope = btn.dataset.scope;
    if (STATE.scope === scope) return;
    STATE.scope = scope;
    document.querySelectorAll('.scope-toggle button').forEach((b) =>
      b.classList.toggle('active', b === btn));
    updateScopeColumnLabels();
    // Re-render the table + re-render the currently-selected variant detail
    renderVariants().then(() => {
      if (STATE.selectedVariant) {
        const tr = document.querySelector(
          `#variants-table tbody tr[data-tag="${STATE.selectedVariant}"]`);
        selectVariant(STATE.selectedVariant, tr);
      }
    });
  })
);

// -------------------- flux sandbox --------------------
let SANDBOX_OVERRIDES = [];
async function initSandbox() {
  const sel = document.getElementById('sandbox-variant');
  if (!sel.dataset.init) {
    const m = await loadManifest();
    sel.innerHTML = '';
    m.variants.forEach((v) => {
      const o = document.createElement('option');
      o.value = v.tag;
      o.textContent = `${v.tag} — ${v.title}`;
      sel.appendChild(o);
    });
    sel.dataset.init = '1';
  }
}

function renderOverrides() {
  const div = document.getElementById('sandbox-overrides');
  if (SANDBOX_OVERRIDES.length === 0) {
    div.innerHTML = '<p class="hint">No overrides — running this will just be the variant\'s reversibility map.</p>';
    return;
  }
  div.innerHTML = SANDBOX_OVERRIDES.map((ov, i) =>
    `<div class="ov-row">
       <code>${escapeHtml(ov.rxn)}</code> → <strong>${escapeHtml(ov.mode)}</strong>
       <button class="small" data-i="${i}">remove</button>
     </div>`).join('');
  div.querySelectorAll('button').forEach((b) =>
    b.addEventListener('click', () => {
      SANDBOX_OVERRIDES.splice(Number(b.dataset.i), 1);
      renderOverrides();
    })
  );
}

document.getElementById('ov-add').addEventListener('click', () => {
  if (STATE.staticMode) return;
  const rxn = document.getElementById('ov-rxn').value.trim();
  const mode = document.getElementById('ov-mode').value;
  if (!rxn) return;
  SANDBOX_OVERRIDES.push({ rxn, mode });
  document.getElementById('ov-rxn').value = '';
  renderOverrides();
});

document.getElementById('sandbox-run').addEventListener('click', async () => {
  if (STATE.staticMode) return;
  const variant = document.getElementById('sandbox-variant').value;
  const modelsRaw = document.getElementById('sandbox-models').value.trim();
  const models = modelsRaw ? modelsRaw.split(/[,\s]+/).filter(Boolean) : null;
  const overrides = {};
  SANDBOX_OVERRIDES.forEach((o) => { overrides[o.rxn] = o.mode; });
  const btn = document.getElementById('sandbox-run');
  const stat = document.getElementById('sandbox-status');
  const out = document.getElementById('sandbox-results');
  btn.disabled = true;
  stat.textContent = 'running FBA on panel…';
  out.innerHTML = '';
  try {
    const t0 = performance.now();
    const res = await API.panelFba({ variant, overrides, models });
    stat.textContent = `done in ${(performance.now() - t0) / 1000 | 0}s (${res.n_models} models, ${res.elapsed_s}s server)`;
    renderSandbox(res, out);
  } catch (exc) {
    stat.textContent = 'error: ' + exc.message;
  } finally {
    btn.disabled = false;
  }
});

function renderSandbox(res, out) {
  const rs = res.results;
  const nGrow = rs.filter((r) => r.grows).length;
  const fluxes = rs.map((r) => r.growth_flux);
  const mean = fluxes.reduce((a, b) => a + b, 0) / Math.max(1, fluxes.length);
  out.innerHTML = `
    <div class="results-summary">
      <div class="badges">
        <div class="badge">variant: <span class="n">${escapeHtml(res.variant)}</span></div>
        <div class="badge">overrides: <span class="n">${res.n_overrides}</span></div>
        <div class="badge">grew: <span class="n">${nGrow}/${rs.length}</span></div>
        <div class="badge">mean growth flux: <span class="n">${mean.toFixed(4)}</span></div>
      </div>
    </div>
    <table class="changed-by-table">
      <thead><tr><th>model_id</th><th>status</th><th>grows</th><th>growth flux</th><th>n overrides applied</th></tr></thead>
      <tbody>${rs.sort((a, b) => b.growth_flux - a.growth_flux).map((r) =>
        `<tr><td>${escapeHtml(r.model_id)}</td>
              <td>${escapeHtml(r.status)}</td>
              <td>${r.grows ? '✓' : '✗'}</td>
              <td class="num">${r.growth_flux.toFixed(4)}</td>
              <td class="num">${r.n_overrides}</td></tr>`
      ).join('')}</tbody>
    </table>`;
}

// -------------------- panel models --------------------
async function loadPanelModelData() {
  if (STATE.panelModels && STATE.panelModelVariants) return;
  const [desc, pmv] = await Promise.all([
    API.data('panel_model_descriptions.json'),
    API.data('panel_model_variants.json'),
  ]);
  STATE.panelModels = desc;
  STATE.panelModelVariants = pmv;
  await loadPanelRxnsets();
  if (!STATE.reactionsPanel) STATE.reactionsPanel = await API.data('reactions_panel.json');
  await loadManifest();
  // Per-reaction perturbation pipeline (precomputed; optional — tolerate absence).
  try { STATE.panelPipeline = await API.data('panel_model_rxn_pipeline.json'); }
  catch (e) { STATE.panelPipeline = null; }
  // Key-reaction direction-sensitivity (precomputed; optional).
  try { STATE.panelKeyReactions = await API.data('panel_key_reactions.json'); }
  catch (e) { STATE.panelKeyReactions = null; }
  // Growth-control knockout essentiality (precomputed; optional).
  try { STATE.panelGrowthControl = await API.data('panel_growth_control.json'); }
  catch (e) { STATE.panelGrowthControl = null; }
  // Synthetic-lethal pairs + flux variability (precomputed; optional).
  try { STATE.panelSyntheticLethal = await API.data('panel_synthetic_lethal.json'); }
  catch (e) { STATE.panelSyntheticLethal = null; }
  try { STATE.panelFva = await API.data('panel_fva.json'); }
  catch (e) { STATE.panelFva = null; }
  // Per-reaction direction options (<,>,=,?) + heuristic calls (precomputed; optional).
  try { STATE.panelRxnDirEffects = await API.data('reaction_direction_effects_panel.json'); }
  catch (e) { STATE.panelRxnDirEffects = null; }
}

function pmEntry(mid, tag) {
  return (STATE.panelModelVariants[mid] || {})[tag] || null;
}
// variant with the largest |Δ flux| for a model -> {tag, e, d} or null
function pmMaxAbs(mid) {
  const m = STATE.panelModelVariants[mid] || {};
  let best = null;
  for (const [tag, e] of Object.entries(m)) {
    const d = Math.abs(e.delta_flux || 0);
    if (!best || d > best.d) best = { tag, e, d };
  }
  return best;
}

async function initPanelModels() {
  const list = document.getElementById('panel-model-list');
  list.innerHTML = '<li><em class="hint">Loading panel models…</em></li>';
  await loadPanelModelData();
  const vsel = document.getElementById('pm-variant');
  if (vsel.options.length <= 1) {
    (await loadManifest()).variants.forEach((v) => {
      if (v.tag === 'baseline') return;
      const o = document.createElement('option');
      o.value = v.tag; o.textContent = `${v.tag} — ${v.title}`;
      vsel.appendChild(o);
    });
  }
  if (!STATE._pmWired) {
    STATE._pmWired = true;
    document.getElementById('pm-search').addEventListener('input', renderPanelModelList);
    document.getElementById('pm-variant').addEventListener('change', (e) => {
      STATE.pmVariantFilter = e.target.value; renderPanelModelList();
    });
    document.getElementById('pm-flux-filter').addEventListener('change', (e) => {
      STATE.pmFluxFilter = e.target.value; renderPanelModelList();
    });
  }
  renderPanelModelList();
}

function renderPanelModelList() {
  const q = (document.getElementById('pm-search').value || '').toLowerCase().trim();
  const varFilt = STATE.pmVariantFilter;
  const fluxFilt = STATE.pmFluxFilter;
  const rows = [];
  for (const mid of Object.keys(STATE.panelModels || {})) {
    const desc = STATE.panelModels[mid] || {};
    let e = null;
    if (varFilt === 'any') { const best = pmMaxAbs(mid); if (best) e = best.e; }
    else e = pmEntry(mid, varFilt);
    const d = e ? (e.delta_flux || 0) : 0;
    const grew = e ? (e.baseline_grows === false && e.variant_grows === true) : false;
    const died = e ? (e.baseline_grows === true && e.variant_grows === false) : false;
    let pass = true;
    if (fluxFilt === 'changed') pass = Math.abs(d) > 1e-6;
    else if (fluxFilt === 'up') pass = d > 1e-6;
    else if (fluxFilt === 'down') pass = d < -1e-6;
    else if (fluxFilt === 'grew') pass = grew;
    else if (fluxFilt === 'died') pass = died;
    if (fluxFilt !== 'any' && !e) pass = false;
    if (!pass) continue;
    if (q && !`${mid} ${desc.organism_name || ''}`.toLowerCase().includes(q)) continue;
    rows.push({ mid, desc, e, d });
  }
  const fluxActive = (varFilt !== 'any') || (fluxFilt !== 'any');
  rows.sort((a, b) => fluxActive ? (Math.abs(b.d) - Math.abs(a.d)) : a.mid.localeCompare(b.mid));
  document.getElementById('pm-result-count').textContent =
    `${rows.length} model${rows.length === 1 ? '' : 's'}` +
    (varFilt !== 'any' ? ` · variant ${varFilt}` : '');
  const ul = document.getElementById('panel-model-list');
  ul.innerHTML = '';
  rows.slice(0, 200).forEach(({ mid, desc, e, d }) => {
    const li = document.createElement('li');
    li.dataset.mid = mid;
    if (mid === STATE.selectedModel) li.classList.add('selected');
    const badge = (e && Math.abs(d) > 1e-6)
      ? ` <span class="pfc-val ${d > 0 ? 'pos' : 'neg'}">${d > 0 ? '+' : ''}${d.toFixed(2)}</span>` : '';
    li.innerHTML = `<strong>${escapeHtml(mid)}</strong>${badge}` +
      `<span class="rxn-name">${escapeHtml(desc.organism_name || '(organism unknown)')}</span>`;
    li.addEventListener('click', () => selectModel(mid));
    ul.appendChild(li);
  });
  if (rows.length > 200) {
    const li = document.createElement('li');
    li.innerHTML = `<em class="hint">… ${rows.length - 200} more (narrow with search)</em>`;
    ul.appendChild(li);
  }
}

function selectModel(mid) {
  STATE.selectedModel = mid;
  STATE.pmOverrides = {};
  document.querySelectorAll('#panel-model-list li').forEach((li) =>
    li.classList.toggle('selected', li.dataset.mid === mid));
  renderPanelModelDetail(mid);
}

function pmTaxChips(desc) {
  const ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'];
  const chips = ranks.map((r) => desc[r]).filter(Boolean)
    .map((v) => `<span class="tax-chip">${escapeHtml(v)}</span>`).join('');
  return chips ? `<div class="tax-chips">${chips}</div>`
               : '<p class="hint">taxonomy unavailable for this assembly</p>';
}

// Diverging Δ-flux bar chart for one model, one row per variant (reuses pfc-* css).
function renderModelVariantChart(rows) {
  const data = rows.filter((e) => e.delta_flux != null);
  const maxAbs = data.length ? Math.max(...data.map((e) => Math.abs(e.delta_flux))) : 0;
  if (!data.length || maxAbs < 1e-6) {
    return `<p class="hint">No variant changes this model's growth flux (all |Δ| &lt; 1e-6).</p>`;
  }
  const sorted = [...data].sort((a, b) => b.delta_flux - a.delta_flux);
  const rowsHtml = sorted.map((e) => {
    const sign = e.delta_flux > 1e-6 ? 'pos' : (e.delta_flux < -1e-6 ? 'neg' : 'zero');
    const half = (Math.abs(e.delta_flux) / maxAbs * 50).toFixed(2);
    const bar = sign === 'zero' ? '' : `<span class="pfc-bar ${sign}" style="width:${half}%"></span>`;
    const flip = (e.baseline_grows != null && e.baseline_grows !== e.variant_grows);
    const icon = flip ? (e.variant_grows ? '<span class="pfc-flip-icon up">↑</span>'
                                         : '<span class="pfc-flip-icon down">↓</span>') : '';
    const tip = `${e.tag}\nΔ growth flux: ${e.delta_flux.toFixed(4)}${flip ? '   ← grow-status flipped' : ''}`;
    return `<div class="pfc-row${flip ? ' flipped' : ''}" title="${escapeHtml(tip)}">` +
      `<span class="pfc-id">${escapeHtml(e.tag)}</span>` +
      `<span class="pfc-flip">${icon}</span>` +
      `<span class="pfc-bar-cell">${bar}</span>` +
      `<span class="pfc-val ${sign}">${e.delta_flux >= 0 ? '+' : ''}${e.delta_flux.toFixed(3)}</span></div>`;
  }).join('');
  return `<div class="panel-flux-chart">
    <div class="pfc-legend">${sorted.length} variants &nbsp;·&nbsp; bar scale ±${maxAbs.toFixed(2)} growth flux units</div>
    <div class="pfc-rows">${rowsHtml}</div></div>`;
}

function renderPanelModelDetail(mid) {
  const pane = document.getElementById('panel-model-detail');
  const desc = STATE.panelModels[mid] || {};
  const nRxn = (STATE.panelRxnsets[mid] || []).length;
  const impactRows = Object.entries(STATE.panelModelVariants[mid] || {})
    .map(([tag, e]) => ({ tag, ...e }))
    .sort((a, b) => Math.abs(b.delta_flux || 0) - Math.abs(a.delta_flux || 0));
  const fmt = (x) => (x == null ? '—' : Number(x).toLocaleString());
  const fmtFlux = (x) => (x == null ? '—' : Number(x).toFixed(3));

  const descCard = `
    <div class="variant-change-card">
      <div class="vcc-label">Organism</div>
      <div class="vcc-title">${escapeHtml(desc.organism_name || '(organism unknown)')}</div>
      ${pmTaxChips(desc)}
      <dl class="pm-desc-dl">
        <dt>reactions</dt><dd>${fmt(desc.n_reactions)}</dd>
        <dt>metabolites</dt><dd>${fmt(desc.n_metabolites)}</dd>
        <dt>genes</dt><dd>${fmt(desc.n_genes)}</dd>
        <dt>open exchanges</dt><dd>${fmt(desc.n_exchanges_open)}</dd>
        <dt>baseline growth flux</dt><dd>${fmtFlux(desc.growth_flux)}</dd>
        <dt>NCBI tax id</dt><dd>${escapeHtml(desc.tax_id || '—')}</dd>
      </dl>
      ${desc.reason ? `<div class="vcc-refs"><span class="vcc-refs-label">in panel because</span><span>${escapeHtml(desc.reason)}</span></div>` : ''}
    </div>`;

  const impactTable = impactRows.length ? `
    <table class="changed-by-table pm-impact-table">
      <thead><tr><th>variant</th><th class="num">Δ growth flux</th><th>grow base→var</th><th class="num">rxns changed in model</th></tr></thead>
      <tbody>${impactRows.map((e) => {
        const sign = e.delta_flux > 1e-6 ? 'pos' : (e.delta_flux < -1e-6 ? 'neg' : 'zero');
        const flip = (e.baseline_grows != null && e.baseline_grows !== e.variant_grows);
        const flag = flip ? (e.variant_grows ? ' ↑' : ' ↓') : '';
        const gb = e.baseline_grows == null ? '—' : (e.baseline_grows ? '✓' : '✗');
        const gv = e.variant_grows == null ? '—' : (e.variant_grows ? '✓' : '✗');
        return `<tr class="pm-impact-row${flip ? ' flipped' : ''}" data-tag="${escapeHtml(e.tag)}">
          <td><span class="tag">${escapeHtml(e.tag)}</span></td>
          <td class="num pfc-val ${sign}">${e.delta_flux >= 0 ? '+' : ''}${Number(e.delta_flux).toFixed(3)}</td>
          <td>${gb} → ${gv}${flag}</td>
          <td class="num">${fmt(e.n_changed)}</td></tr>`;
      }).join('')}</tbody>
    </table>
    <p class="hint">Click a variant row to list the reactions it changes in this model.</p>
    <div id="pm-variant-rxns"></div>`
    : '<p class="hint">No variant changes any reaction present in this model.</p>';

  const liveHtml = `
    <div class="live-only">
      <h3>Reaction directions &amp; manual override <span class="hint">— live FBA on ${escapeHtml(mid)}</span></h3>
      <p class="hint">Every unique reaction in this model: its <strong>default</strong> (baseline) direction,
        the selected variant's direction (changes highlighted), and a dropdown to set your own. Then run FBA
        on this model to see the growth effect of your overrides vs baseline.</p>
      <div class="sandbox-controls">
        <label>Compare to variant: <select id="pm-ov-variant"></select></label>
        <button id="pm-ov-apply" class="small">set overrides = this variant's changes</button>
      </div>
      <div id="pm-overrides" class="overrides-box"></div>
      <div class="search-row">
        <input id="pm-ov-search" type="text" placeholder="filter this model's reactions…" autocomplete="off">
        <span id="pm-ov-rxn-count" class="hint"></span>
      </div>
      <div id="pm-ov-rxnlist" class="pm-ov-rxnlist"></div>
      <div class="run-row">
        <button id="pm-ov-run" class="primary">Run FBA on this model</button>
        <span id="pm-ov-status" class="hint"></span>
      </div>
      <div id="pm-ov-results"></div>
    </div>
    <p class="hint static-only">Manual per-model FBA requires live mode (<code>python3 site/serve.py --live</code>).</p>`;

  pane.innerHTML = `
    <h3>${escapeHtml(mid)} <span class="hint">— ${escapeHtml(desc.organism_name || '(organism unknown)')}</span></h3>
    ${descCard}
    <h3>Reactions &amp; impact</h3>
    <div class="flux-impact-grid">
      <div class="card"><h4>unique reactions</h4><div class="stat">${nRxn.toLocaleString()}</div></div>
      <div class="card"><h4>variants changing ≥1 rxn here</h4><div class="stat">${impactRows.filter((e) => e.n_changed > 0).length}</div></div>
      <div class="card"><h4>variants that move growth flux</h4><div class="stat">${impactRows.filter((e) => Math.abs(e.delta_flux) > 1e-6).length}</div></div>
    </div>
    <h3>Per-variant impact on this model <span class="hint">— Δ growth flux vs baseline; click a row for the changed reactions</span></h3>
    ${impactTable}
    <h3>Δ growth flux by variant <span class="hint">— this model, each variant vs the baseline cascade</span></h3>
    ${renderModelVariantChart(impactRows)}
    <h3>Reaction-direction heuristics — growth under &lt;, &gt;, =, ?
      <span class="hint">— starting from the default, each reaction is set one-at-a-time to each option; where the 4 heuristics send it, and the resulting growth (? = knocked off)</span></h3>
    <div class="rde-controls">
      <input id="pm-rde-search" type="text" placeholder="filter reactions…" autocomplete="off">
      <label><input type="checkbox" id="pm-rde-sensitive" checked> only growth-sensitive</label>
      <label><input type="checkbox" id="pm-rde-disagree"> only heuristic disagreements</label>
      <span id="pm-rde-count" class="hint"></span>
    </div>
    <div id="pm-rde-charts"></div>
    ${liveHtml}
    <h3>Heuristic perturbation pipeline <span class="hint">— single-reaction &amp; cumulative growth-flux effect of one heuristic, vs baseline</span></h3>
    <div class="pm-pipe-controls">
      <label>Heuristic: <select id="pm-pipe-variant"></select></label>
      <span id="pm-pipe-note" class="hint"></span>
    </div>
    <div id="pm-pipe-charts"></div>
    <h3>Key reactions — direction sensitivity <span class="hint">— reactions whose direction change most disrupts this model's growth (every reaction probed in all directions vs baseline)</span></h3>
    <div id="pm-key-charts"></div>
    <h3>Growth control — knockout essentiality <span class="hint">— blocking each reaction (FBA reaction deletion): which reactions keep growth high (essential) or low (limiting)</span></h3>
    <div id="pm-gc-charts"></div>
    <h3>Synthetic-lethal pairs <span class="hint">— reaction pairs whose joint knockout collapses growth though neither alone does</span></h3>
    <div id="pm-sl-charts"></div>
    <h3>Flux variability <span class="hint">— at near-optimal growth: which reactions are blocked, flux-forced (obligate), or flexible</span></h3>
    <div id="pm-fva-charts"></div>`;

  pane.querySelectorAll('.pm-impact-row').forEach((tr) =>
    tr.addEventListener('click', () => {
      pane.querySelectorAll('.pm-impact-row').forEach((r) => r.classList.remove('selected'));
      tr.classList.add('selected');
      showModelVariantRxns(mid, tr.dataset.tag);
    }));

  if (!STATE.staticMode) initPmOverride(mid);
  initPmPipeline(mid);        // precomputed (static) data — render in both modes
  renderPmKeyReactions(mid);   // precomputed direction-sensitivity — both modes
  renderPmGrowthControl(mid);  // precomputed knockout essentiality + limiting metabolites
  renderPmSyntheticLethal(mid);// precomputed synthetic-lethal pairs
  renderPmFva(mid);            // precomputed flux variability
  renderPmRxnDirEffects(mid);  // per-reaction direction options (<,>,=,?) + 4 heuristic calls
}

// ----- key reactions: per-reaction direction sensitivity (precomputed) -----
// Ranks reactions by how much flipping their direction (vs baseline) moves growth.
// Data: site/data/panel_key_reactions.json from scripts/build_key_reactions.py.
function renderPmKeyReactions(mid) {
  const host = document.getElementById('pm-key-charts');
  if (!host) return;
  const data = STATE.panelKeyReactions;
  const m = (data && data.models) ? data.models[mid] : null;
  if (!m) {
    host.innerHTML = `<p class="hint">Key-reaction data not available (run <code>scripts/build_key_reactions.py</code>).</p>`;
    return;
  }
  const rxns = m.reactions || [];
  if (!rxns.length) {
    host.innerHTML = `<p class="hint">No single reaction's direction change moves this model's growth.</p>`;
    return;
  }
  host.innerHTML = `
    <p class="hint">baseline growth flux ${Number(m.base_flux).toFixed(3)} &nbsp;·&nbsp;
      ${m.n_tested} reactions probed in every direction &nbsp;·&nbsp;
      ${rxns.length} shown (|Δ growth| &gt; 1e-6, ranked by severity)</p>
    ${renderKeyReactionChart(rxns)}
    ${renderGlobalKeyCard((data && data.global) || [])}`;
  bindRxnLinks(host);
}

// Diverging severity bar chart for one model's key reactions (reuses pfc-* styles).
function renderKeyReactionChart(rxns) {
  const maxSev = Math.max(...rxns.map((r) => r.severity), 1e-9);
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = rxns.map((r) => {
    const d = r.best_delta;
    const sign = d > 1e-6 ? 'pos' : (d < -1e-6 ? 'neg' : 'zero');
    const half = (r.severity / maxSev * 50).toFixed(2);
    const bar = `<span class="pfc-bar ${sign}" style="width:${half}%"></span>`;
    const byd = Object.entries(r.by_dir)
      .map(([k, v]) => `${k} ${v >= 0 ? '+' : ''}${Number(v).toFixed(2)}`).join('   ');
    const tip = `${r.rxn}${rxnName(r.rxn) ? ' — ' + rxnName(r.rxn) : ''}\n` +
      `baseline ${r.base_dir} → ${r.best_dir}\nΔ growth by forced direction: ${byd}`;
    return `<div class="pfc-row" title="${escapeHtml(tip)}">` +
      `<span class="pfc-id"><a href="#" class="rxn-link" data-rxn="${escapeHtml(r.rxn)}">${escapeHtml(r.rxn)}</a></span>` +
      `<span class="pfc-keydir">${revBadge(r.base_dir)}→${revBadge(r.best_dir)}</span>` +
      `<span class="pfc-bar-cell">${bar}</span>` +
      `<span class="pfc-val ${sign}">${d >= 0 ? '+' : ''}${Number(d).toFixed(3)}</span></div>`;
  }).join('');
  return `<div class="panel-flux-chart key-rxn-chart">
    <div class="pfc-legend">${rxns.length} key reactions &nbsp;·&nbsp; bar = severity (max |Δ growth| over the forced directions) &nbsp;·&nbsp;
      <span class="pfc-val neg">red</span> = a direction change kills growth, <span class="pfc-val pos">green</span> = boosts it</div>
    <div class="pfc-rows">${rows}</div></div>`;
}

// Collapsible cross-panel tally: reactions that are key in the most models.
function renderGlobalKeyCard(glob) {
  if (!glob || !glob.length) return '';
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = glob.slice(0, 15).map((g) =>
    `<tr><td><a href="#" class="rxn-link" data-rxn="${escapeHtml(g.rxn)}">${escapeHtml(g.rxn)}</a></td>` +
    `<td>${escapeHtml(rxnName(g.rxn))}</td>` +
    `<td class="num">${g.n_models}</td>` +
    `<td class="num">${Number(g.max_severity).toFixed(2)}</td>` +
    `<td class="num">${Number(g.mean_severity).toFixed(2)}</td></tr>`).join('');
  return `<details class="key-rxn-global">
    <summary>Most frequently key across the 100 panel models</summary>
    <table class="changed-by-table">
      <thead><tr><th>rxn</th><th>name</th><th># models</th><th>max |Δ|</th><th>mean |Δ|</th></tr></thead>
      <tbody>${rows}</tbody></table></details>`;
}

// ----- reaction-direction heuristics: growth under <,>,=,? per reaction -----
// Data: site/data/reaction_direction_effects_panel.json (build_reaction_direction_effects.py).
// For each reaction (set one-at-a-time from the model's default bounds) we show the
// growth under each of the 4 options and where the 4 heuristics send it.
const RDE_OPTS = [
  { k: '<', lab: '&lt; reverse' }, { k: '>', lab: '&gt; forward' },
  { k: '=', lab: '= reversible' }, { k: '?', lab: '? off' },
];
const RDE_SCHEMES = [
  { k: 'default', i: 'D', lab: 'Default (model)' }, { k: 'jankowski', i: 'J', lab: 'Jankowski (group contribution)' },
  { k: 'flamholz', i: 'F', lab: 'Flamholz 2012 (eQuilibrator)' }, { k: 'opus', i: 'O', lab: 'Claude Opus 4.8' },
];
function rdeColor(g, base) {
  if (g == null) return '#eee';
  const r = base > 1e-6 ? g / base : (g > 1e-6 ? 1 : 0);
  if (r < 0.02) return '#f7c7c2';
  if (r < 0.5) return '#f6d9a8';
  if (r < 0.98) return '#f2eeb0';
  if (r <= 1.02) return '#cbe7c8';
  return '#b9d4ef';
}
function rdeDirHtml(d) { return d === '<' ? '&lt;' : d === '>' ? '&gt;' : d; }
function rdeDirClass(d) {
  return d === '<' ? 'lt' : d === '>' ? 'gt' : d === '=' ? 'eq' : d === '?' ? 'q' : 'na';
}
function rdeSensitive(rec) {
  const v = Object.values(rec.g).filter((x) => x != null);
  return v.length ? (Math.max(...v) - Math.min(...v)) > 1e-6 : false;
}
function rdeDisagree(rec) {
  const calls = RDE_SCHEMES.map((s) => rec.dirs[s.k]).filter((c) => c && c !== 'NA' && c !== '?');
  return new Set(calls).size > 1;
}
function renderPmRxnDirEffects(mid) {
  const host = document.getElementById('pm-rde-charts');
  if (!host) return;
  const data = STATE.panelRxnDirEffects;
  const m = (data && data.models) ? data.models[mid] : null;
  if (!m) {
    host.innerHTML = '<p class="hint">Reaction-direction-effects data not available '
      + '(run <code>scripts/build_reaction_direction_effects.py</code>).</p>';
    return;
  }
  STATE._rde = { mid, base: m.base_flux, rows: m.reactions };
  renderRdeTable();
  ['pm-rde-search', 'pm-rde-sensitive', 'pm-rde-disagree'].forEach((id) => {
    const el = document.getElementById(id);
    if (el) { el.addEventListener('input', renderRdeTable); el.addEventListener('change', renderRdeTable); }
  });
}
function renderRdeTable() {
  const host = document.getElementById('pm-rde-charts');
  const st = STATE._rde;
  if (!host || !st) return;
  const q = (document.getElementById('pm-rde-search').value || '').trim().toLowerCase();
  const onlySen = document.getElementById('pm-rde-sensitive').checked;
  const onlyDis = document.getElementById('pm-rde-disagree').checked;
  const rxnName = (id) => (STATE.reactionsPanel && STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  let rows = st.rows;
  if (onlySen) rows = rows.filter(rdeSensitive);
  if (onlyDis) rows = rows.filter(rdeDisagree);
  if (q) rows = rows.filter((r) => r.rxn.toLowerCase().includes(q) || rxnName(r.rxn).toLowerCase().includes(q));
  const shown = rows.slice(0, 300);
  const ob = STATE.panelRxnDirEffects.option_bounds || {};
  const head = '<tr><th>reaction</th><th>heuristic calls</th>'
    + RDE_OPTS.map((o) => `<th class="rde-opt" title="bounds ${JSON.stringify(ob[o.k] || '')}">${o.lab}</th>`).join('')
    + '</tr>';
  const body = shown.map((rec) => {
    const badges = RDE_SCHEMES.map((s) => {
      const d = rec.dirs[s.k];
      return `<span class="rde-badge rde-d-${rdeDirClass(d)}" title="${s.lab}: ${d}">${s.i}:${rdeDirHtml(d)}</span>`;
    }).join(' ');
    const cells = RDE_OPTS.map((o) => {
      const g = rec.g[o.k];
      const who = RDE_SCHEMES.filter((s) => rec.dirs[s.k] === o.k)
        .map((s) => `<span class="rde-who" title="${s.lab} sends it here">${s.i}</span>`).join('');
      const val = g == null ? '×' : (Number.isInteger(g) ? g : g.toFixed(g < 10 ? 2 : 1));
      return `<td class="rde-cell" style="background:${rdeColor(g, st.base)}">`
        + `<span class="rde-val">${val}</span>${who ? `<span class="rde-whos">${who}</span>` : ''}</td>`;
    }).join('');
    return `<tr><td class="rde-rxn"><a href="#" class="rxn-link" data-rxn="${escapeHtml(rec.rxn)}">${escapeHtml(rec.rxn)}</a>`
      + `${rxnName(rec.rxn) ? `<span class="rde-name">${escapeHtml(rxnName(rec.rxn))}</span>` : ''}</td>`
      + `<td class="rde-heur">${badges}</td>${cells}</tr>`;
  }).join('');
  const cnt = document.getElementById('pm-rde-count');
  if (cnt) cnt.textContent = `${rows.length} reactions${rows.length > shown.length ? ` (showing ${shown.length})` : ''}`
    + ` · baseline growth ${Number(st.base).toFixed(3)}`;
  host.innerHTML =
    '<div class="rde-legend">Each cell = biomass growth when that one reaction is forced to the option (all others left at default). '
    + '<b style="background:#cbe7c8">≈baseline</b> <b style="background:#b9d4ef">boosted</b> '
    + '<b style="background:#f2eeb0">reduced</b> <b style="background:#f6d9a8">strongly ↓</b> '
    + '<b style="background:#f7c7c2">~dead</b> <b style="background:#eee">× infeasible</b>. '
    + 'Letters in a cell = which heuristics send the reaction there (D default, J Jankowski, F Flamholz, O Opus). '
    + '<b>?</b> = unknown, tested as off/knockout.</div>'
    + `<div class="rde-tablewrap"><table class="rde-table"><thead>${head}</thead><tbody>${body}</tbody></table></div>`
    + renderRdeGlobal();
  bindRxnLinks(host);
}
function renderRdeGlobal() {
  const g = STATE.panelRxnDirEffects && STATE.panelRxnDirEffects.global;
  if (!g) return '';
  const row = (x) => `<tr><td><a href="#" class="rxn-link" data-rxn="${x.base}">${x.base}</a></td><td class="num">${x.n_models}</td></tr>`;
  const dis = (g.most_disagreed || []).slice(0, 12).map(row).join('');
  const off = (g.most_off_essential || []).slice(0, 12).map(row).join('');
  return '<details class="rde-global"><summary>Cross-panel patterns (all 100 models)</summary>'
    + '<div class="rde-global-grid">'
    + `<div><h5>Heuristics disagree most</h5><table class="changed-by-table"><thead><tr><th>rxn</th><th># models</th></tr></thead><tbody>${dis}</tbody></table></div>`
    + `<div><h5>Essential when off (? = 0 growth)</h5><table class="changed-by-table"><thead><tr><th>rxn</th><th># models</th></tr></thead><tbody>${off}</tbody></table></div>`
    + '</div></details>';
}

async function showModelVariantRxns(mid, tag) {
  const out = document.getElementById('pm-variant-rxns');
  out.innerHTML = '<p class="loading">loading…</p>';
  const p = await loadVariant(tag);
  const modelRxns = new Set(STATE.panelRxnsets[mid] || []);
  const changed = (p.diffs || []).filter((d) => modelRxns.has(d.rxn))
    .sort((a, b) => a.rxn.localeCompare(b.rxn));
  const e = pmEntry(mid, tag) || {};
  if (!changed.length) {
    out.innerHTML = `<p class="hint">Variant <span class="tag">${escapeHtml(tag)}</span> changes no reaction present in this model.</p>`;
    return;
  }
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const dsign = e.delta_flux > 1e-6 ? 'pos' : (e.delta_flux < -1e-6 ? 'neg' : 'zero');
  out.innerHTML = `
    <div class="pm-variant-rxns-head">
      <span class="tag">${escapeHtml(tag)}</span> changes <strong>${changed.length}</strong>
      reaction${changed.length === 1 ? '' : 's'} in this model &nbsp;·&nbsp; model Δ flux
      <span class="pfc-val ${dsign}">${(e.delta_flux || 0) >= 0 ? '+' : ''}${Number(e.delta_flux || 0).toFixed(3)}</span>
    </div>
    <table class="changed-by-table">
      <thead><tr><th>rxn</th><th>baseline</th><th>variant</th><th>name</th></tr></thead>
      <tbody>${changed.map((d) =>
        `<tr><td><a href="#" class="rxn-link" data-rxn="${escapeHtml(d.rxn)}">${escapeHtml(d.rxn)}</a></td>
             <td>${revBadge(d.base)}</td><td>${revBadge(d.new)}</td>
             <td>${escapeHtml(rxnName(d.rxn))}</td></tr>`).join('')}</tbody>
    </table>`;
  bindRxnLinks(out);
}

// ----- per-reaction perturbation pipeline (precomputed; both modes) -----
// For one model + one heuristic: each variant-changed reaction flipped alone vs
// baseline (marginal), then applied cumulatively in |Δ| order. Data precomputed
// in site/data/panel_model_rxn_pipeline.json by scripts/build_panel_rxn_pipeline.py.
function initPmPipeline(mid) {
  const sel = document.getElementById('pm-pipe-variant');
  const charts = document.getElementById('pm-pipe-charts');
  if (!sel || !charts) return;
  const data = (STATE.panelPipeline || {})[mid] || {};
  const tags = Object.keys(data);
  if (!tags.length) {
    sel.innerHTML = '';
    sel.disabled = true;
    charts.innerHTML = STATE.panelPipeline
      ? `<p class="hint">No heuristic changes the direction of any reaction in this model, so there is nothing to perturb.</p>`
      : `<p class="hint">Perturbation pipeline data not available (run <code>scripts/build_panel_rxn_pipeline.py</code>).</p>`;
    const note = document.getElementById('pm-pipe-note');
    if (note) note.textContent = '';
    return;
  }
  // Order the dropdown by |full Δ| (cumulative endpoint) descending.
  const fullDelta = (t) => {
    const c = data[t].cumulative;
    return c.length ? Math.abs(c[c.length - 1].delta) : 0;
  };
  tags.sort((a, b) => fullDelta(b) - fullDelta(a));
  const titles = (STATE.manifest && STATE.manifest.variants || [])
    .reduce((m, v) => { m[v.tag] = v.title; return m; }, {});
  sel.disabled = false;
  sel.innerHTML = tags.map((t) =>
    `<option value="${escapeHtml(t)}">${escapeHtml(t)}${titles[t] ? ' — ' + escapeHtml(titles[t]) : ''}</option>`
  ).join('');
  // Default: the top-filter variant if it has data here, else the biggest mover.
  let def = (STATE.pmVariantFilter && STATE.pmVariantFilter !== 'any' && data[STATE.pmVariantFilter])
    ? STATE.pmVariantFilter : tags[0];
  sel.value = def;
  sel.onchange = () => renderPmPipeline(mid, sel.value);
  renderPmPipeline(mid, def);
}

function renderPmPipeline(mid, tag) {
  const charts = document.getElementById('pm-pipe-charts');
  if (!charts) return;
  const d = ((STATE.panelPipeline || {})[mid] || {})[tag];
  if (!d) { charts.innerHTML = `<p class="hint">No perturbation data for this variant.</p>`; return; }
  const note = document.getElementById('pm-pipe-note');
  if (note) {
    note.textContent =
      `baseline growth flux ${Number(d.base_flux).toFixed(3)} · ` +
      `${d.singles.length} reaction${d.singles.length === 1 ? '' : 's'} changed in this model`;
  }
  charts.innerHTML = `
    <h4>Single-reaction Δ growth flux <span class="hint">— each changed reaction applied alone vs baseline, sorted by |Δ|</span></h4>
    ${renderPmMarginalChart(d.singles)}
    <h4>Cumulative Δ growth flux <span class="hint">— the ranked changes applied 1, then 1+2, … vs baseline</span></h4>
    ${renderPmCumulativeChart(d.cumulative)}`;
  bindRxnLinks(charts);
}

// Sorted diverging bar chart of single-reaction Δ growth flux (reuses pfc-* styles).
function renderPmMarginalChart(singles) {
  if (!singles || !singles.length) return '<p class="hint">No changed reactions.</p>';
  const maxAbs = Math.max(...singles.map((s) => Math.abs(s.delta)), 1e-9);
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = singles.map((s) => {
    const sign = s.delta > 1e-6 ? 'pos' : (s.delta < -1e-6 ? 'neg' : 'zero');
    const halfPct = (Math.abs(s.delta) / maxAbs * 50).toFixed(2);
    const bar = sign === 'zero' ? '' : `<span class="pfc-bar ${sign}" style="width:${halfPct}%"></span>`;
    const tip = `${s.rxn}${rxnName(s.rxn) ? ' — ' + rxnName(s.rxn) : ''}\n` +
      `single-change Δ growth flux: ${s.delta >= 0 ? '+' : ''}${Number(s.delta).toFixed(3)}`;
    return `<div class="pfc-row" title="${escapeHtml(tip)}">` +
      `<span class="pfc-id"><a href="#" class="rxn-link" data-rxn="${escapeHtml(s.rxn)}">${escapeHtml(s.rxn)}</a></span>` +
      `<span class="pfc-flip"></span>` +
      `<span class="pfc-bar-cell">${bar}</span>` +
      `<span class="pfc-val ${sign}">${s.delta >= 0 ? '+' : ''}${Number(s.delta).toFixed(3)}</span></div>`;
  }).join('');
  const nNonzero = singles.filter((s) => Math.abs(s.delta) > 1e-6).length;
  return `<div class="panel-flux-chart">
    <div class="pfc-legend">${singles.length} changed reaction${singles.length === 1 ? '' : 's'} &nbsp;·&nbsp;
      ${nNonzero} individually move growth flux &nbsp;·&nbsp; bar scale ±${maxAbs.toFixed(2)} growth flux units</div>
    <div class="pfc-rows">${rows}</div></div>`;
}

// Inline SVG line chart of cumulative Δ growth flux as ranked changes are applied.
function renderPmCumulativeChart(cumulative) {
  const n = (cumulative || []).length;
  if (!n) return '<p class="hint">No changed reactions.</p>';
  const W = 600, H = 210, padL = 60, padR = 18, padT = 14, padB = 30;
  const xs = (i) => padL + (n === 1 ? (W - padL - padR) / 2 : (i / (n - 1)) * (W - padL - padR));
  const deltas = cumulative.map((c) => c.delta);
  let lo = Math.min(0, ...deltas), hi = Math.max(0, ...deltas);
  if (hi === lo) { hi += 1; lo -= 1; }
  const ys = (v) => padT + (1 - (v - lo) / (hi - lo)) * (H - padT - padB);
  const y0 = ys(0);
  const pts = cumulative.map((c, i) => `${xs(i).toFixed(1)},${ys(c.delta).toFixed(1)}`).join(' ');
  const dots = cumulative.map((c, i) => {
    const sign = c.delta > 1e-6 ? 'pos' : (c.delta < -1e-6 ? 'neg' : 'zero');
    const tip = `after ${i + 1} reaction${i ? 's' : ''} (＋${c.rxn}): Δ ${c.delta >= 0 ? '+' : ''}${Number(c.delta).toFixed(3)}`;
    return `<circle class="pm-cum-dot ${sign}" cx="${xs(i).toFixed(1)}" cy="${ys(c.delta).toFixed(1)}" r="2.6"><title>${escapeHtml(tip)}</title></circle>`;
  }).join('');
  const last = cumulative[n - 1];
  const lsign = last.delta > 1e-6 ? 'pos' : (last.delta < -1e-6 ? 'neg' : 'zero');
  return `<div class="pm-cum-chart">
    <svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="cumulative delta growth flux">
      <line class="pm-cum-axis" x1="${padL}" y1="${y0.toFixed(1)}" x2="${(W - padR).toFixed(1)}" y2="${y0.toFixed(1)}"/>
      <text class="pm-cum-lbl" x="${padL - 7}" y="${(ys(hi) + 4).toFixed(1)}" text-anchor="end">${hi.toFixed(1)}</text>
      <text class="pm-cum-lbl" x="${padL - 7}" y="${(y0 + 4).toFixed(1)}" text-anchor="end">0</text>
      <text class="pm-cum-lbl" x="${padL - 7}" y="${(ys(lo) + 4).toFixed(1)}" text-anchor="end">${lo.toFixed(1)}</text>
      <polyline class="pm-cum-line" points="${pts}"/>
      ${dots}
      <text class="pm-cum-lbl" x="${padL}" y="${H - 9}" text-anchor="start">1</text>
      <text class="pm-cum-lbl" x="${(W - padR).toFixed(1)}" y="${H - 9}" text-anchor="end">${n}</text>
      <text class="pm-cum-axtitle" x="${((padL + W - padR) / 2).toFixed(1)}" y="${H - 9}" text-anchor="middle"># reactions applied (in |Δ| rank order)</text>
    </svg>
    <div class="pfc-legend">cumulative Δ growth flux vs baseline as the top-ranked changes are applied; endpoint = full variant Δ =
      <span class="pfc-val ${lsign}">${last.delta >= 0 ? '+' : ''}${Number(last.delta).toFixed(3)}</span></div>
  </div>`;
}

// ----- growth control: knockout essentiality + flux at optimum (precomputed) -----
function renderPmGrowthControl(mid) {
  const host = document.getElementById('pm-gc-charts');
  if (!host) return;
  const gc = STATE.panelGrowthControl;
  const m = (gc && gc.models) ? gc.models[mid] : null;
  if (!m) {
    host.innerHTML = `<p class="hint">Growth-control data not available (run <code>scripts/build_growth_control.py</code>).</p>`;
    return;
  }
  const rx = m.reactions || [];
  if (!rx.length) {
    host.innerHTML = `<p class="hint">No single reaction knockout changes this model's growth.</p>`;
    return;
  }
  host.innerHTML = `
    <p class="hint">baseline growth flux ${Number(m.base_flux).toFixed(3)} &nbsp;·&nbsp; ${m.n_tested} reactions knocked out &nbsp;·&nbsp;
      <span class="pfc-val neg">${m.n_essential}</span> essential (KO collapses growth) &nbsp;·&nbsp;
      <span class="pfc-val pos">${m.n_limiting}</span> growth-limiting (KO raises growth)</p>
    <h4>Knockout impact <span class="hint">— Δ growth when each reaction is blocked, top by magnitude</span></h4>
    ${renderGcBars(rx.slice(0, 40))}
    <h4>Flux vs essentiality <span class="hint">— flux carried at the growth optimum vs knockout Δ growth (each reaction)</span></h4>
    ${renderGcScatter(rx)}
    <h4>Limiting metabolites <span class="hint">— LP shadow prices: metabolite pools that most constrain growth</span></h4>
    ${renderLimitingMetaboliteBars(m.metabolites || [])}`;
  bindRxnLinks(host);
}

function renderLimitingMetaboliteBars(mets) {
  if (!mets.length) return '<p class="hint">No non-zero metabolite shadow prices.</p>';
  const top = mets.slice(0, 18);
  const maxAbs = Math.max(...top.map((x) => Math.abs(x.shadow_price)), 1e-9);
  const rows = top.map((x) => {
    const sp = x.shadow_price, sign = sp < 0 ? 'neg' : 'pos';
    const w = (Math.abs(sp) / maxAbs * 100).toFixed(1);
    return `<div class="an-distrow" title="${escapeHtml(x.name)} (${x.met}): shadow price ${sp >= 0 ? '+' : ''}${Number(sp).toFixed(4)}">` +
      `<span class="an-distlbl">${escapeHtml(x.name)}</span>` +
      `<span class="an-distbar"><span class="an-seg" style="width:${w}%;background:var(--accent-2)"></span></span>` +
      `<span class="an-distn ${sign}">${sp >= 0 ? '+' : ''}${Number(sp).toFixed(3)}</span></div>`;
  }).join('');
  return `<div class="an-distwrap an-crit">${rows}</div>
    <p class="hint">|shadow price| = how strongly that metabolite's mass balance limits growth at the optimum.</p>`;
}

// Synthetic-lethal pairs: ranked table of double-knockouts that collapse growth.
function renderPmSyntheticLethal(mid) {
  const host = document.getElementById('pm-sl-charts');
  if (!host) return;
  const sl = STATE.panelSyntheticLethal;
  const m = (sl && sl.models) ? sl.models[mid] : null;
  if (!m) { host.innerHTML = `<p class="hint">Synthetic-lethal data not available (run <code>scripts/build_synthetic_lethal.py</code>).</p>`; return; }
  if (!m.pairs || !m.pairs.length) {
    host.innerHTML = `<p class="hint">No synthetic-lethal pairs found among the top ${m.n_candidates} flux-carrying non-essential reactions in this model.</p>`;
    return;
  }
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = m.pairs.map((p) =>
    `<tr><td><a href="#" class="rxn-link" data-rxn="${escapeHtml(p.a)}">${escapeHtml(p.a)}</a> + <a href="#" class="rxn-link" data-rxn="${escapeHtml(p.b)}">${escapeHtml(p.b)}</a></td>` +
    `<td class="num pfc-val neg">${Number(p.joint_delta).toFixed(2)}</td>` +
    `<td class="num">${Number(p.single_a).toFixed(2)} / ${Number(p.single_b).toFixed(2)}</td>` +
    `<td class="num pfc-val ${p.epistasis < -1e-6 ? 'neg' : 'zero'}">${Number(p.epistasis).toFixed(2)}</td>` +
    `<td title="${escapeHtml(rxnName(p.a) + ' + ' + rxnName(p.b))}">${escapeHtml((rxnName(p.a) || '?').slice(0, 22))} + …</td></tr>`).join('');
  host.innerHTML = `<p class="hint">${m.n_pairs} synthetic-lethal/sick pair${m.n_pairs === 1 ? '' : 's'} among the top ${m.n_candidates} flux-carrying, individually-non-essential reactions.</p>
    <table class="changed-by-table">
      <thead><tr><th>reaction pair</th><th>joint Δ growth</th><th>singles (A / B)</th><th>epistasis</th><th>names</th></tr></thead>
      <tbody>${rows}</tbody></table>
    <p class="hint">Joint knockout collapses growth though each single barely moves it; <strong>epistasis</strong> = joint − (single A + single B), strongly negative = synergistic lethality.</p>`;
  bindRxnLinks(host);
}

// Flux variability: classification summary + flux-range bars for the key reactions.
function renderPmFva(mid) {
  const host = document.getElementById('pm-fva-charts');
  if (!host) return;
  const fv = STATE.panelFva;
  const m = (fv && fv.models) ? fv.models[mid] : null;
  if (!m) { host.innerHTML = `<p class="hint">Flux-variability data not available (run <code>scripts/build_fva.py</code>).</p>`; return; }
  const rx = (m.reactions || []).filter((r) => r.kind !== 'blocked').slice(0, 36);
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  // shared symmetric scale across shown reactions
  const ext = Math.max(...rx.map((r) => Math.max(Math.abs(r.min), Math.abs(r.max))), 1e-9);
  const W = 720, rowH = 16, padL = 130, padR = 50, top = 8;
  const H = top + rx.length * rowH + 8;
  const X = (v) => padL + ((v + ext) / (2 * ext)) * (W - padL - padR);
  const zero = X(0);
  const bars = rx.map((r, i) => {
    const y = top + i * rowH;
    const x1 = X(r.min), x2 = X(r.max);
    const col = r.kind === 'flux_forced' ? _rgba(_C_WARN, 0.85) : _rgba(_C_ACC, 0.7);
    return `<a href="#" class="rxn-link" data-rxn="${escapeHtml(r.rxn)}"><text class="an-lbl" x="${padL - 6}" y="${y + rowH - 5}" text-anchor="end">${escapeHtml(r.rxn)}</text></a>` +
      `<rect x="${Math.min(x1, x2).toFixed(1)}" y="${y + 2}" width="${Math.max(2, Math.abs(x2 - x1)).toFixed(1)}" height="${rowH - 6}" fill="${col}"><title>${escapeHtml(r.rxn + (rxnName(r.rxn) ? ' — ' + rxnName(r.rxn) : ''))}: [${Number(r.min).toFixed(1)}, ${Number(r.max).toFixed(1)}] (${r.kind})</title></rect>`;
  }).join('');
  host.innerHTML = `<p class="hint">baseline growth flux ${Number(m.base_flux).toFixed(3)} &nbsp;·&nbsp;
    <span class="pfc-val neg">${m.n_forced}</span> flux-forced (obligate) &nbsp;·&nbsp;
    ${m.n_flexible} flexible &nbsp;·&nbsp; <span class="hint">${m.n_blocked}</span> blocked (cannot carry flux)</p>
    <div class="an-scroll an-tall"><svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">
      <line class="an-grid" x1="${zero.toFixed(1)}" y1="2" x2="${zero.toFixed(1)}" y2="${H - 4}"/>
      <text class="an-lbl" x="${padL}" y="${H - 2}" text-anchor="start">${(-ext).toFixed(0)}</text>
      <text class="an-lbl" x="${(W - padR).toFixed(1)}" y="${H - 2}" text-anchor="end">+${ext.toFixed(0)}</text>
      ${bars}</svg></div>
    <p class="hint">Flux interval each reaction can take at ≥99% of optimal growth (Mahadevan &amp; Schilling 2003). <span class="pfc-val neg">red</span> = flux-forced (interval excludes 0 → obligate for growth); teal = flexible. Blocked reactions omitted.</p>`;
  bindRxnLinks(host);
}

function renderGcBars(rx) {
  const maxAbs = Math.max(...rx.map((r) => Math.abs(r.ko_delta)), 1e-9);
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = rx.map((r) => {
    const d = r.ko_delta;
    const sign = d < -1e-6 ? 'neg' : (d > 1e-6 ? 'pos' : 'zero');
    const half = (Math.abs(d) / maxAbs * 50).toFixed(2);
    const bar = `<span class="pfc-bar ${sign}" style="width:${half}%"></span>`;
    const tip = `${r.rxn}${rxnName(r.rxn) ? ' — ' + rxnName(r.rxn) : ''}\n` +
      `knockout Δ growth: ${d >= 0 ? '+' : ''}${Number(d).toFixed(3)} (${r.kind})\n` +
      `flux at optimum: ${Number(r.flux_opt).toFixed(2)}   reduced cost: ${Number(r.reduced_cost).toFixed(3)}`;
    return `<div class="pfc-row" title="${escapeHtml(tip)}">` +
      `<span class="pfc-id"><a href="#" class="rxn-link" data-rxn="${escapeHtml(r.rxn)}">${escapeHtml(r.rxn)}</a></span>` +
      `<span class="pfc-keydir">${r.kind === 'essential' ? '⛔' : r.kind === 'limiting' ? '▲' : '·'}</span>` +
      `<span class="pfc-bar-cell">${bar}</span>` +
      `<span class="pfc-val ${sign}">${d >= 0 ? '+' : ''}${Number(d).toFixed(3)}</span></div>`;
  }).join('');
  return `<div class="panel-flux-chart key-rxn-chart">
    <div class="pfc-legend"><span class="pfc-val neg">red ⛔</span> = essential (blocking collapses growth) &nbsp;·&nbsp; <span class="pfc-val pos">green ▲</span> = growth-limiting (blocking raises growth) &nbsp;·&nbsp; bar = |Δ growth|</div>
    <div class="pfc-rows">${rows}</div></div>`;
}

function renderGcScatter(rx) {
  const W = 720, H = 380, padL = 64, padB = 46, padT = 16, padR = 18;
  const fmax = Math.max(...rx.map((r) => Math.abs(r.flux_opt)), 1);
  const ys = rx.map((r) => r.ko_delta);
  let lo = Math.min(0, ...ys), hi = Math.max(0, ...ys);
  if (hi <= lo) hi = lo + 1;
  const X = (f) => padL + Math.sqrt(f / fmax) * (W - padL - padR);  // sqrt: flux spans 0..~1000
  const Y = (v) => padT + (1 - (v - lo) / (hi - lo)) * (H - padT - padB);
  const col = (k) => k === 'essential' ? _rgba(_C_WARN, 0.75) : k === 'limiting' ? _rgba(_C_GOOD, 0.75) : _rgba(_C_ACC, 0.4);
  const rxnName = (id) => (STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const y0 = Y(0);
  const dots = rx.map((r) => {
    const tip = `${r.rxn}${rxnName(r.rxn) ? ' — ' + rxnName(r.rxn) : ''}: flux ${Number(r.flux_opt).toFixed(1)}, KO Δ ${r.ko_delta >= 0 ? '+' : ''}${Number(r.ko_delta).toFixed(2)} (${r.kind})`;
    return `<circle cx="${X(Math.abs(r.flux_opt)).toFixed(1)}" cy="${Y(r.ko_delta).toFixed(1)}" r="3.5" fill="${col(r.kind)}" stroke="rgba(0,0,0,0.25)" stroke-width="0.5"><title>${escapeHtml(tip)}</title></circle>`;
  }).join('');
  const xt = [0, 0.25, 0.5, 1].map((f) => {
    const fv = f * fmax;
    return `<text class="an-lbl" x="${X(fv).toFixed(1)}" y="${H - padB + 14}" text-anchor="middle">${Math.round(fv)}</text>`;
  }).join('');
  return `<div class="an-scroll"><svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">
    <line class="an-grid" x1="${padL}" y1="${y0.toFixed(1)}" x2="${W - padR}" y2="${y0.toFixed(1)}"/>
    <text class="an-lbl" x="${padL - 6}" y="${(Y(hi) + 3).toFixed(1)}" text-anchor="end">${hi.toFixed(1)}</text>
    <text class="an-lbl" x="${padL - 6}" y="${(y0 + 3).toFixed(1)}" text-anchor="end">0</text>
    <text class="an-lbl" x="${padL - 6}" y="${(Y(lo) + 3).toFixed(1)}" text-anchor="end">${lo.toFixed(1)}</text>
    ${xt}${dots}
    <text class="an-axt" x="${(padL + W - padR) / 2}" y="${H - 4}" text-anchor="middle">|flux| carried at growth optimum (√ scale)</text>
    <text class="an-axt" transform="translate(15,${(padT + H - padB) / 2}) rotate(-90)" text-anchor="middle">knockout Δ growth</text></svg></div>
    <p class="hint">Bottom-left/red (very negative Δ) = essential bottlenecks; bottom-right/red = high-flux essential backbone; on the zero line = dispensable; above zero/green = growth-limiting.</p>`;
}

// ----- per-model live override panel (live mode only) -----
async function initPmOverride(mid) {
  const vsel = document.getElementById('pm-ov-variant');
  if (vsel) {
    vsel.innerHTML = '';
    (await loadManifest()).variants.forEach((v) => {
      if (v.tag === 'baseline') return;
      const o = document.createElement('option');
      o.value = v.tag; o.textContent = `${v.tag} — ${v.title}`;
      vsel.appendChild(o);
    });
    // Default the comparison to the variant selected in the top filter (if any),
    // else the variant with the biggest flux impact on this model.
    let def = (STATE.pmVariantFilter && STATE.pmVariantFilter !== 'any') ? STATE.pmVariantFilter : null;
    if (!def) { const best = pmMaxAbs(mid); def = best ? best.tag : (vsel.options[0] && vsel.options[0].value); }
    if (def) vsel.value = def;
  }
  if (!STATE.baselineMap) {
    try { STATE.baselineMap = (await API.data('baseline.json')).map; }
    catch (e) { STATE.baselineMap = {}; }
  }
  STATE.pmOverrides = {};
  await loadCompareVariant(vsel ? vsel.value : null);
  renderPmOverrides();
  renderPmOvRxnList(mid);
  if (vsel) vsel.onchange = async () => {
    await loadCompareVariant(vsel.value);
    renderPmOvRxnList(STATE.selectedModel);
  };
  const search = document.getElementById('pm-ov-search');
  if (search) search.oninput = () => renderPmOvRxnList(STATE.selectedModel);
  const runBtn = document.getElementById('pm-ov-run');
  if (runBtn) runBtn.onclick = () => runPmOverride(mid);
  const applyBtn = document.getElementById('pm-ov-apply');
  if (applyBtn) applyBtn.onclick = () => applyVariantOverrides(mid);
}

// Load the {rxn: new_dir} diff map for the variant being compared against baseline.
async function loadCompareVariant(tag) {
  STATE.pmCompareTag = tag;
  STATE.pmCompareDiffs = {};
  if (!tag || tag === 'baseline') return;
  const p = await loadVariant(tag);
  for (const d of (p.diffs || [])) STATE.pmCompareDiffs[d.rxn] = d.new;
}

// Pre-fill the override dropdowns with the compared variant's direction for every
// reaction it changes in this model (so you can run that heuristic, then tweak).
function applyVariantOverrides(mid) {
  const REV_TO_MODE = { '>': 'forward', '<': 'reverse', '=': 'free' };
  const base = STATE.baselineMap || {};
  const inModel = new Set(STATE.panelRxnsets[mid] || []);
  for (const [rxn, nv] of Object.entries(STATE.pmCompareDiffs || {})) {
    if (inModel.has(rxn) && nv !== base[rxn] && REV_TO_MODE[nv]) {
      STATE.pmOverrides[rxn] = REV_TO_MODE[nv];
    }
  }
  renderPmOverrides();
  renderPmOvRxnList(mid);
}

function renderPmOverrides() {
  const div = document.getElementById('pm-overrides');
  if (!div) return;
  const keys = Object.keys(STATE.pmOverrides || {});
  if (!keys.length) {
    div.innerHTML = '<p class="hint">No overrides set — choose reaction directions below.</p>';
    return;
  }
  div.innerHTML = keys.map((rxn) =>
    `<div class="ov-row"><code>${escapeHtml(rxn)}</code> → <strong>${escapeHtml(STATE.pmOverrides[rxn])}</strong>
       <button class="small" data-rxn="${escapeHtml(rxn)}">remove</button></div>`).join('');
  div.querySelectorAll('button').forEach((b) =>
    b.addEventListener('click', () => {
      delete STATE.pmOverrides[b.dataset.rxn];
      renderPmOverrides();
      renderPmOvRxnList(STATE.selectedModel);
    }));
}

function renderPmOvRxnList(mid) {
  const wrap = document.getElementById('pm-ov-rxnlist');
  if (!wrap) return;
  const q = (document.getElementById('pm-ov-search').value || '').toLowerCase().trim();
  let rxns = STATE.panelRxnsets[mid] || [];
  if (q) rxns = rxns.filter((id) =>
    id.toLowerCase().includes(q) ||
    ((STATE.reactionsPanel[id]?.name || '').toLowerCase().includes(q)));
  const base = STATE.baselineMap || {};
  const diffs = STATE.pmCompareDiffs || {};
  const tag = STATE.pmCompareTag;
  const varDir = (id) => (id in diffs ? diffs[id] : base[id]);
  const nChanged = rxns.filter((id) => varDir(id) !== base[id]).length;
  const count = document.getElementById('pm-ov-rxn-count');
  if (count) count.textContent =
    `${rxns.length} reactions${tag ? ` · ${nChanged} changed by ${tag}` : ''}`;
  const modes = ['as_is', 'forward', 'reverse', 'free', 'off'];
  const varHdr = tag ? escapeHtml(tag) : 'variant';
  const rowsHtml = rxns.slice(0, 200).map((id) => {
    const d = base[id] || '?';
    const v = varDir(id) || d;
    const changed = v !== d;
    const sel = STATE.pmOverrides[id] || 'as_is';
    const name = STATE.reactionsPanel[id]?.name || '';
    return `<tr class="pm-ov-rxn-row${changed ? ' pm-changed' : ''}${sel !== 'as_is' ? ' active' : ''}">
      <td><a href="#" class="rxn-link" data-rxn="${escapeHtml(id)}">${escapeHtml(id)}</a></td>
      <td class="pm-ov-name" title="${escapeHtml(name)}">${escapeHtml(name)}</td>
      <td>${revBadge(d)}</td>
      <td>${revBadge(v)}${changed ? ' <span class="pm-chg" title="changed by this variant">●</span>' : ''}</td>
      <td><select data-rxn="${escapeHtml(id)}">${modes.map((mm) =>
        `<option value="${mm}"${mm === sel ? ' selected' : ''}>${mm}</option>`).join('')}</select></td>
    </tr>`;
  }).join('');
  wrap.innerHTML = `<table class="changed-by-table pm-ov-table">
      <thead><tr><th>reaction</th><th>name</th><th>default</th><th>${varHdr}</th><th>override</th></tr></thead>
      <tbody>${rowsHtml}</tbody></table>` +
    (rxns.length > 200 ? `<p class="hint">… ${rxns.length - 200} more (filter to narrow)</p>` : '');
  wrap.querySelectorAll('select').forEach((s) =>
    s.addEventListener('change', () => {
      const rxn = s.dataset.rxn;
      if (s.value === 'as_is') delete STATE.pmOverrides[rxn];
      else STATE.pmOverrides[rxn] = s.value;
      renderPmOverrides();
      s.closest('.pm-ov-rxn-row').classList.toggle('active', s.value !== 'as_is');
    }));
  bindRxnLinks(wrap);
}

async function runPmOverride(mid) {
  const overrides = { ...STATE.pmOverrides };
  const btn = document.getElementById('pm-ov-run');
  const status = document.getElementById('pm-ov-status');
  const out = document.getElementById('pm-ov-results');
  if (!Object.keys(overrides).length) {
    out.innerHTML = '<p class="hint">Set at least one reaction’s override (or use '
      + '“set overrides = this variant’s changes”) before running.</p>';
    return;
  }
  btn.disabled = true; status.textContent = 'running FBA…'; out.innerHTML = '';
  try {
    const t0 = performance.now();
    // Always referenced to the baseline cascade so it matches the "default"
    // column; the variant selector only drives the comparison column above.
    const [base, mod] = await Promise.all([
      API.panelFba({ variant: 'baseline', overrides: {}, models: [mid], n_workers: 1 }),
      API.panelFba({ variant: 'baseline', overrides, models: [mid], n_workers: 1 }),
    ]);
    const b = (base.results || [])[0] || {};
    const v = (mod.results || [])[0] || {};
    const d = (v.growth_flux || 0) - (b.growth_flux || 0);
    const sign = d > 1e-6 ? 'pos' : (d < -1e-6 ? 'neg' : 'zero');
    status.textContent = `done in ${((performance.now() - t0) / 1000) | 0}s`;
    out.innerHTML = `
      <div class="flux-impact-grid">
        <div class="card"><h4>baseline</h4>
          <div class="stat">${(b.growth_flux || 0).toFixed(4)}</div>
          <div class="stat-sub">${b.grows ? 'grows' : 'no growth'}</div></div>
        <div class="card"><h4>with ${Object.keys(overrides).length} override(s)</h4>
          <div class="stat">${(v.growth_flux || 0).toFixed(4)}</div>
          <div class="stat-sub">${v.grows ? 'grows' : 'no growth'}${v.status ? ` · ${escapeHtml(v.status)}` : ''}</div></div>
        <div class="card"><h4>Δ growth flux</h4>
          <div class="stat pfc-val ${sign}">${d >= 0 ? '+' : ''}${d.toFixed(4)}</div>
          <div class="stat-sub">${b.grows !== v.grows ? (v.grows ? 'became grower ↑' : 'stopped growing ↓') : 'grow-status unchanged'}</div></div>
      </div>`;
  } catch (exc) {
    status.textContent = 'error: ' + exc.message;
  } finally { btn.disabled = false; }
}

// -------------------- resizable / collapsible list panes --------------------
// Insert a drag handle between each .reaction-grid's list pane and detail pane so
// the sidebar can be resized with the cursor (and double-click to collapse/expand).
// Width persists across tabs/reloads via localStorage.
function enhanceResizableGrids() {
  const KEY = 'reaction-grid-list-w';
  const saved = localStorage.getItem(KEY);
  document.querySelectorAll('.reaction-grid').forEach((grid) => {
    if (grid.dataset.resizable) return;
    grid.dataset.resizable = '1';
    const listPane = grid.children[0];
    const detailPane = grid.children[1];
    if (!listPane || !detailPane) return;
    if (saved) grid.style.setProperty('--list-w', saved);

    const resizer = document.createElement('div');
    resizer.className = 'grid-resizer';
    resizer.title = 'Drag to resize · double-click to collapse/expand';
    resizer.innerHTML = '<span class="grid-resizer-grip"></span>';
    grid.insertBefore(resizer, detailPane);

    const curW = () => getComputedStyle(grid).getPropertyValue('--list-w').trim() || '360px';
    let dragging = false, startX = 0, startW = 0;
    const onMove = (e) => {
      if (!dragging) return;
      const w = Math.max(0, Math.min(startW + (e.clientX - startX), grid.clientWidth - 200));
      grid.style.setProperty('--list-w', w + 'px');
      grid.classList.toggle('list-collapsed', w < 8);
    };
    const onUp = () => {
      if (!dragging) return;
      dragging = false;
      document.body.style.userSelect = '';
      localStorage.setItem(KEY, curW());
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseup', onUp);
    };
    resizer.addEventListener('mousedown', (e) => {
      dragging = true;
      startX = e.clientX;
      startW = listPane.getBoundingClientRect().width;
      document.body.style.userSelect = 'none';
      window.addEventListener('mousemove', onMove);
      window.addEventListener('mouseup', onUp);
      e.preventDefault();
    });
    resizer.addEventListener('dblclick', () => {
      const nowCollapsed = grid.classList.toggle('list-collapsed');
      if (nowCollapsed) {
        const w = curW();
        grid._prevListW = (w && w !== '0px') ? w : '360px';
        grid.style.setProperty('--list-w', '0px');
      } else {
        grid.style.setProperty('--list-w', grid._prevListW || '360px');
      }
      localStorage.setItem(KEY, curW());
    });
  });
}

// -------------------- analytics tab --------------------
const _C_GOOD = [106, 209, 122], _C_WARN = [239, 111, 108], _C_ACC = [95, 216, 184];
const _rgba = (c, a) => `rgba(${c[0]},${c[1]},${c[2]},${Math.max(0, Math.min(1, a)).toFixed(3)})`;
const _seq = (t) => _rgba(_C_ACC, Math.max(0, Math.min(1, t)) * 0.9 + 0.07);

async function initAnalytics() {
  STATE.analyticsInit = true;
  try { STATE.methodCmp = await API.data('method_comparison.json'); } catch (e) { STATE.methodCmp = null; }
  if (!STATE.panelKeyReactions) { try { STATE.panelKeyReactions = await API.data('panel_key_reactions.json'); } catch (e) {} }
  if (!STATE.panelGrowthControl) { try { STATE.panelGrowthControl = await API.data('panel_growth_control.json'); } catch (e) {} }
  if (!STATE.panelSyntheticLethal) { try { STATE.panelSyntheticLethal = await API.data('panel_synthetic_lethal.json'); } catch (e) {} }
  if (!STATE.panelFva) { try { STATE.panelFva = await API.data('panel_fva.json'); } catch (e) {} }
  if (!STATE.reactionsPanel) { try { STATE.reactionsPanel = await API.data('reactions_panel.json'); } catch (e) {} }

  // --- direction-method comparison (KEGG_default reference vs heuristics; 4 reaction scopes) ---
  // scope toggle: re-render only the method matrices (panel analyses are scope-independent)
  document.querySelectorAll('#an-scope-toggle button').forEach((btn) =>
    btn.addEventListener('click', () => {
      const s = btn.dataset.mscope;
      if (STATE.methodScope === s) return;
      STATE.methodScope = s;
      document.querySelectorAll('#an-scope-toggle button').forEach((b) =>
        b.classList.toggle('active', b === btn));
      renderMethodMatrices();
    }));
  renderMethodMatrices();

  // --- growth-control panel analyses (100-model panel; scope-independent) ---
  renderKeyCriticality(document.getElementById('an-criticality'));
  renderEssentialityGlobal(document.getElementById('an-essential'));
  renderFvaGlobal(document.getElementById('an-fva'));
  renderSLGlobal(document.getElementById('an-sl'));
  renderLimitingMetabolitesGlobal(document.getElementById('an-metabolites'));
}

// --- comparison matrices (agreement, confusion, direction makeup), N methods ---
// All read method_comparison.json for the currently-selected reaction scope; the
// method set is per-scope (KEGG_default + heuristics in model scopes, heuristics-only wide).
function renderMethodMatrices() {
  const mc = STATE.methodCmp;
  const agHost = document.getElementById('an-agreement');
  const cfHost = document.getElementById('an-confusion');
  const ddHost = document.getElementById('an-dirdist');
  if (!mc || !mc.modes) {
    [agHost, cfHost, ddHost].forEach((h) => { if (h) h.innerHTML = '<p class="hint">Method-comparison data unavailable (run build_method_comparison.py).</p>'; });
    return;
  }
  const mode = mc.modes[STATE.methodScope] || mc.modes.all;
  const nEl = document.getElementById('an-mc-n');
  if (nEl) {
    const n = (mc.counts && mc.counts[STATE.methodScope]) || mode.n_reactions;
    nEl.textContent = Number(n).toLocaleString();
  }
  renderMethodAgreement(mode, agHost);
  renderMethodConfusion(mode, cfHost);
  renderMethodDirDist(mode, ddHost);
}

// Direction symbols → human labels + colors, shared by the method matrices.
const _DIR_LABEL = { '>': 'forward', '<': 'reverse', '=': 'reversible', '?': 'unknown' };
const _DIR_COLOR = { '>': 'var(--good)', '<': 'var(--warn)', '=': 'var(--accent-2)', '?': 'var(--border)' };

// 1. Method × method directional-agreement heatmap (3×3, rank-similarity ordered).
function renderMethodAgreement(mode, host) {
  if (!host) return;
  const methods = mode.methods, n = methods.length, A = mode.agreement, N = mode.agreement_n;
  // rescale color over the off-diagonal range so the few cells read distinctly
  let amin = 1;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i !== j) amin = Math.min(amin, A[i][j]);
  const resc = (a) => (1 - amin > 1e-6 ? (a - amin) / (1 - amin) : 1);
  const cs = 76, padL = 132, padT = 118, W = padL + n * cs + 8, H = padT + n * cs + 8;
  let cells = '', rowlbl = '', collbl = '';
  for (let i = 0; i < n; i++) {
    const y = padT + i * cs;
    rowlbl += `<text class="an-lbl" x="${padL - 8}" y="${y + cs / 2 + 3}" text-anchor="end">${escapeHtml(methods[i])}</text>`;
    collbl += `<text class="an-lbl" transform="translate(${padL + i * cs + cs / 2},${padT - 8}) rotate(-30)" text-anchor="start">${escapeHtml(methods[i])}</text>`;
    for (let j = 0; j < n; j++) {
      const a = A[i][j], x = padL + j * cs, t = resc(a);
      const tip = i === j
        ? `${methods[i]}: makes a directional call (>, <, =) on ${N[i][i].toLocaleString()} reactions`
        : `${methods[i]} vs ${methods[j]}: ${(a * 100).toFixed(1)}% call the same direction over ${N[i][j].toLocaleString()} reactions both decide`;
      cells += `<rect x="${x}" y="${y}" width="${cs - 2}" height="${cs - 2}" rx="3" fill="${_rgba(_C_ACC, i === j ? 1 : 0.14 + 0.84 * t)}"><title>${escapeHtml(tip)}</title></rect>`;
      const dark = (i === j || t > 0.5);
      cells += `<text x="${x + cs / 2}" y="${y + cs / 2 - 2}" text-anchor="middle" font-family="var(--mono)" font-size="17" font-weight="600" fill="${dark ? '#10141a' : '#cdd6e2'}" pointer-events="none">${Math.round(a * 100)}%</text>`;
      cells += `<text x="${x + cs / 2}" y="${y + cs / 2 + 14}" text-anchor="middle" font-family="var(--mono)" font-size="9" fill="${dark ? '#2a3340' : '#9aa7b8'}" pointer-events="none">n=${N[i][j].toLocaleString()}</text>`;
    }
  }
  host.innerHTML = `<div class="an-scroll"><svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">${collbl}${rowlbl}${cells}</svg></div>
    <p class="hint">Fraction of reactions <strong>both</strong> methods call directionally (&gt;, &lt;, =) that get the same call; shared "unknown" (?) is excluded. Teal intensity rescaled ${(amin * 100).toFixed(0)}–100% to spread the cells. Methods ordered by rank similarity — adjacent = most alike.</p>`;
}

// 2. Pairwise direction × direction confusion matrices (categories rank-ordered).
function renderMethodConfusion(mode, host) {
  if (!host) return;
  const conf = mode.confusion || [];
  if (!conf.length) { host.innerHTML = '<p class="hint">No confusion data.</p>'; return; }
  const grids = conf.map((c) => {
    // Drop any direction category whose row AND column are entirely empty
    // (e.g. methods that never emit "?") so the grid shows no phantom rows.
    const rowSum = (i) => c.matrix[i].reduce((s, v) => s + v, 0);
    const colSum = (j) => c.matrix.reduce((s, row) => s + row[j], 0);
    const keep = c.cats.map((_, i) => i).filter((i) => rowSum(i) > 0 || colSum(i) > 0);
    const cats = keep.map((i) => c.cats[i]);
    const m = keep.map((i) => keep.map((j) => c.matrix[i][j]));
    const k = cats.length;
    let cmax = 1;
    for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) cmax = Math.max(cmax, m[i][j]);
    // Directional agreement consistent with the agreement heatmap: matches over
    // reactions BOTH methods call directionally (>, <, =), excluding any "?".
    const dir = cats.map((cat, i) => (cat !== '?' ? i : -1)).filter((i) => i >= 0);
    let bothDecide = 0, sameDir = 0;
    for (const i of dir) { for (const j of dir) bothDecide += m[i][j]; sameDir += m[i][i]; }
    const dirAgree = bothDecide ? sameDir / bothDecide : 0;
    const oneSidedUnknown = c.n - bothDecide; // co-covered reactions where one side is "?"
    const cs = 38, padL = 40, padT = 24, W = padL + k * cs + 6, H = padT + k * cs + 6;
    let cells = '', rowlbl = '', collbl = '';
    for (let i = 0; i < k; i++) {
      const y = padT + i * cs;
      rowlbl += `<text class="an-lbl" x="${padL - 6}" y="${y + cs / 2 + 3}" text-anchor="end">${escapeHtml(cats[i])}</text>`;
      collbl += `<text class="an-lbl" x="${padL + i * cs + cs / 2}" y="${padT - 6}" text-anchor="middle">${escapeHtml(cats[i])}</text>`;
      for (let j = 0; j < k; j++) {
        const v = m[i][j], x = padL + j * cs, t = v / cmax;
        const tip = `${c.a} = ${cats[i]} (${_DIR_LABEL[cats[i]]}), ${c.b} = ${cats[j]} (${_DIR_LABEL[cats[j]]}): ${v.toLocaleString()} reactions`;
        cells += `<rect x="${x}" y="${y}" width="${cs - 1.5}" height="${cs - 1.5}" rx="2" fill="${i === j ? _rgba(_C_GOOD, 0.18 + 0.74 * t) : _seq(t)}"><title>${escapeHtml(tip)}</title></rect>`;
        if (v) cells += `<text class="an-cell" x="${x + cs / 2}" y="${y + cs / 2 + 3}" text-anchor="middle" fill="${t > 0.5 ? '#10141a' : '#cdd6e2'}">${v >= 10000 ? Math.round(v / 1000) + 'k' : v}</text>`;
      }
    }
    const sub = `${(dirAgree * 100).toFixed(1)}% same direction · n=${bothDecide.toLocaleString()} both decide`
      + (oneSidedUnknown > 0 ? ` · ${oneSidedUnknown.toLocaleString()} one-sided "?"` : '');
    return `<div class="an-confcell">
      <div class="an-conftitle">${escapeHtml(c.a)} <span class="hint">(rows)</span> × ${escapeHtml(c.b)} <span class="hint">(cols)</span></div>
      <div class="an-confsub">${sub}</div>
      <svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">${collbl}${rowlbl}${cells}</svg></div>`;
  }).join('');
  const legend = `<span class="an-legitem"><span class="an-sw" style="background:${_rgba(_C_GOOD, 0.85)}"></span>same call (diagonal)</span>`
    + `<span class="an-legitem"><span class="an-sw" style="background:${_rgba(_C_ACC, 0.85)}"></span>different call — darker = more reactions</span>`;
  host.innerHTML = `<div class="an-legend">${legend}</div><div class="an-confgrid">${grids}</div>
    <p class="hint">Each cell = # reactions where the row method made the row call (&gt; forward, &lt; reverse, = reversible, ? unknown) and the column method made the column call. The <span class="pfc-val pos">green diagonal</span> = same call; brighter off-diagonal (teal) = larger systematic disagreement (e.g. one method calls reversible "=" where another commits to a direction). The "% same direction" matches the agreement heatmap (over reactions both decide); reactions where one method says "?" are listed separately. Counts ≥ 10,000 abbreviated "k" — exact value on hover. Direction categories are ordered to put the most-confused pair adjacent.</p>`;
}

// 3. Per-method direction makeup: 100% stacked bars (forward / reverse / reversible / unknown).
function renderMethodDirDist(mode, host) {
  if (!host) return;
  const methods = mode.methods, dist = mode.dist;
  const order = ['>', '<', '=', '?'];
  const rows = methods.map((mth) => {
    const d = dist[mth] || {};
    const tot = order.reduce((s, k) => s + (d[k] || 0), 0) || 1;
    const bars = order.map((k) => {
      const v = d[k] || 0, pct = v / tot * 100;
      if (pct < 0.01) return '';
      return `<span class="an-seg" style="width:${pct.toFixed(2)}%;background:${_DIR_COLOR[k]}" title="${escapeHtml(mth)} ${_DIR_LABEL[k]} (${k}): ${v.toLocaleString()} (${pct.toFixed(1)}%)"></span>`;
    }).join('');
    return `<div class="an-distrow"><span class="an-distlbl" title="${escapeHtml(mth)}">${escapeHtml(mth)}</span><span class="an-distbar">${bars}</span><span class="an-distn">${tot.toLocaleString()}</span></div>`;
  }).join('');
  const legend = order.map((k) =>
    `<span class="an-legitem"><span class="an-sw" style="background:${_DIR_COLOR[k]}"></span>${_DIR_LABEL[k]} (${k})</span>`).join('');
  host.innerHTML = `<div class="an-legend">${legend}</div><div class="an-distwrap">${rows}</div>
    <p class="hint">Composition of each method's direction calls, each normalized to 100% over <strong>its own coverage</strong> — the count at right shows how many reactions that method calls, and the method sets only partly overlap (Flamholz covers far fewer reactions than Opus), so read this as each method's internal makeup, not a like-for-like per-reaction comparison. <strong>KEGG_default</strong> (Core models scope) is the direction the models were built with — its split of forward / reverse / reversible is the bar the heuristics are trying to match. Opus 4.8 commits to a forward direction far more often; the ΔG′ methods leave many reactions reversible "=".</p>`;
}

// 10. Essential reactions across the panel (knockout): frequency bars + per-model count histogram.
function renderEssentialityGlobal(host) {
  if (!host) return;
  const gc = STATE.panelGrowthControl;
  if (!gc || !gc.global || !gc.global.length) { host.innerHTML = '<p class="hint">Growth-control data unavailable (run build_growth_control.py).</p>'; return; }
  const rxnName = (id) => (STATE.reactionsPanel && STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const top = gc.global.slice(0, 20);
  const nmax = Math.max(...top.map((x) => x.n_essential_models), 1);
  const bars = top.map((x) =>
    `<div class="an-distrow"><span class="an-distlbl"><a href="#" class="rxn-link" data-rxn="${escapeHtml(x.rxn)}">${escapeHtml(x.rxn)}</a></span>` +
    `<span class="an-distbar"><span class="an-seg" style="width:${(x.n_essential_models / nmax * 100).toFixed(1)}%;background:var(--warn)" title="${escapeHtml(x.rxn)}${rxnName(x.rxn) ? ' — ' + rxnName(x.rxn) : ''}: essential in ${x.n_essential_models} models, max |KO Δ| ${x.max_abs_ko}"></span></span>` +
    `<span class="an-distn">${x.n_essential_models}</span></div>`).join('');
  // histogram: # essential reactions per model
  const counts = Object.values(gc.models).map((m) => m.n_essential || 0);
  const hmax = Math.max(...counts, 1), nb = 20, bw = Math.ceil((hmax + 1) / nb);
  const bins = new Array(nb).fill(0);
  counts.forEach((c) => { bins[Math.min(nb - 1, Math.floor(c / bw))]++; });
  const bmax = Math.max(...bins, 1);
  const W = 460, H = 240, padL = 40, padB = 34, padT = 12, padR = 10;
  const barw = (W - padL - padR) / nb;
  const hbars = bins.map((c, i) => {
    const bh = (c / bmax) * (H - padT - padB);
    return `<rect x="${(padL + i * barw).toFixed(1)}" y="${(H - padB - bh).toFixed(1)}" width="${(barw - 1).toFixed(1)}" height="${bh.toFixed(1)}" fill="${_rgba(_C_WARN, 0.8)}"><title>${i * bw}–${(i + 1) * bw - 1} essential rxns: ${c} models</title></rect>`;
  }).join('');
  host.innerHTML = `<div class="an-twocol">
    <div><h4>Essential in the most models (top 20)</h4><div class="an-distwrap an-crit">${bars}</div></div>
    <div><h4># essential reactions per model</h4><div class="an-scroll"><svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">
      <text class="an-lbl" x="${padL - 4}" y="${padT + 6}" text-anchor="end">${bmax}</text>
      <text class="an-lbl" x="${padL}" y="${H - 6}" text-anchor="start">0</text>
      <text class="an-lbl" x="${W - padR}" y="${H - 6}" text-anchor="end">${hmax}</text>
      ${hbars}
      <text class="an-axt" x="${(padL + W - padR) / 2}" y="${H - 4}" text-anchor="middle"># essential reactions in a model</text></svg></div></div>
  </div>`;
  bindRxnLinks(host);
}

// 6. Reaction criticality from the key-reaction sweep: frequency bars + scatter.
function renderKeyCriticality(host) {
  if (!host) return;
  const kr = STATE.panelKeyReactions;
  if (!kr || !kr.global || !kr.global.length) { host.innerHTML = '<p class="hint">Key-reaction data unavailable (run build_key_reactions.py).</p>'; return; }
  const g = kr.global;
  const rxnName = (id) => (STATE.reactionsPanel && STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const top = g.slice(0, 20);
  const nmax = Math.max(...top.map((x) => x.n_models), 1);
  const bars = top.map((x) => {
    const pct = (x.n_models / nmax * 100).toFixed(1);
    return `<div class="an-distrow"><span class="an-distlbl"><a href="#" class="rxn-link" data-rxn="${escapeHtml(x.rxn)}">${escapeHtml(x.rxn)}</a></span>` +
      `<span class="an-distbar"><span class="an-seg" style="width:${pct}%;background:var(--warn)" title="${escapeHtml(x.rxn)}${rxnName(x.rxn) ? ' — ' + rxnName(x.rxn) : ''}: key in ${x.n_models} models, max |Δ| ${x.max_severity}"></span></span>` +
      `<span class="an-distn">${x.n_models}</span></div>`;
  }).join('');
  // scatter: x=n_models, y=max_severity
  const W = 720, H = 360, padL = 56, padB = 44, padT = 14, padR = 20;
  const xmax = Math.max(...g.map((x) => x.n_models), 1), ymax = Math.max(...g.map((x) => x.max_severity), 1);
  const X = (v) => padL + (v / xmax) * (W - padL - padR);
  const Y = (v) => H - padB - (v / ymax) * (H - padB - padT);
  const dots = g.map((x) =>
    `<circle class="an-bub" cx="${X(x.n_models).toFixed(1)}" cy="${Y(x.max_severity).toFixed(1)}" r="4"><title>${escapeHtml(x.rxn)}${rxnName(x.rxn) ? ' — ' + rxnName(x.rxn) : ''}: key in ${x.n_models} models, max |Δ growth| ${x.max_severity}, mean ${x.mean_severity}</title></circle>`).join('');
  const grid = [0, 0.5, 1].map((f) =>
    `<text class="an-lbl" x="${padL - 6}" y="${(Y(f * ymax) + 3).toFixed(1)}" text-anchor="end">${(f * ymax).toFixed(0)}</text>` +
    `<text class="an-lbl" x="${X(f * xmax).toFixed(1)}" y="${H - padB + 14}" text-anchor="middle">${Math.round(f * xmax)}</text>`).join('');
  host.innerHTML = `<div class="an-twocol">
    <div><h4>Most frequently key (top 20)</h4><div class="an-distwrap an-crit">${bars}</div></div>
    <div><h4>Frequency vs severity</h4><div class="an-scroll"><svg viewBox="0 0 ${W} ${H}" class="an-svg" style="max-width:${W}px">${grid}${dots}
      <text class="an-axt" x="${(padL + W - padR) / 2}" y="${H - 4}" text-anchor="middle"># panel models where reaction is key</text>
      <text class="an-axt" transform="translate(14,${(H - padB) / 2}) rotate(-90)" text-anchor="middle">max |Δ growth|</text></svg></div></div>
  </div>`;
  bindRxnLinks(host);
}

// 11. Limiting metabolites across the panel (shadow prices).
function renderLimitingMetabolitesGlobal(host) {
  if (!host) return;
  const mg = STATE.panelGrowthControl && STATE.panelGrowthControl.metabolites_global;
  if (!mg || !mg.length) { host.innerHTML = '<p class="hint">Shadow-price data unavailable.</p>'; return; }
  const top = mg.slice(0, 20);
  const nmax = Math.max(...top.map((x) => x.n_models), 1);
  const rows = top.map((x) =>
    `<div class="an-distrow" title="${escapeHtml(x.name)} (${x.met}): limiting in ${x.n_models} models, max |shadow price| ${x.max_abs_sp}">` +
    `<span class="an-distlbl">${escapeHtml(x.name)}</span>` +
    `<span class="an-distbar"><span class="an-seg" style="width:${(x.n_models / nmax * 100).toFixed(1)}%;background:var(--accent-2)"></span></span>` +
    `<span class="an-distn">${x.n_models}</span></div>`).join('');
  host.innerHTML = `<div class="an-distwrap an-crit">${rows}</div>
    <p class="hint"># panel models in which each metabolite is among the top growth-limiting pools (largest |shadow price|).</p>`;
}

// 12. Synthetic-lethal pairs recurring across the panel.
function renderSLGlobal(host) {
  if (!host) return;
  const g = STATE.panelSyntheticLethal && STATE.panelSyntheticLethal.global;
  if (!g || !g.length) { host.innerHTML = '<p class="hint">No recurring synthetic-lethal pairs.</p>'; return; }
  const rxnName = (id) => (STATE.reactionsPanel && STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const rows = g.slice(0, 18).map((x) =>
    `<tr><td><a href="#" class="rxn-link" data-rxn="${escapeHtml(x.a)}">${escapeHtml(x.a)}</a> + <a href="#" class="rxn-link" data-rxn="${escapeHtml(x.b)}">${escapeHtml(x.b)}</a></td>` +
    `<td class="num">${x.n_models}</td><td class="num pfc-val neg">${Number(x.max_abs_joint).toFixed(1)}</td>` +
    `<td title="${escapeHtml(rxnName(x.a) + ' + ' + rxnName(x.b))}">${escapeHtml((rxnName(x.a) || '?').slice(0, 26))}…</td></tr>`).join('');
  host.innerHTML = `<table class="changed-by-table">
    <thead><tr><th>reaction pair</th><th># models</th><th>max |joint Δ|</th><th>name (A)</th></tr></thead>
    <tbody>${rows}</tbody></table>
    <p class="hint">Pairs that are synthetic-lethal/sick in the most panel models — conserved jointly-essential reaction couples.</p>`;
  bindRxnLinks(host);
}

// 13. Flux-forced / blocked reactions across the panel (FVA).
function renderFvaGlobal(host) {
  if (!host) return;
  const g = STATE.panelFva && STATE.panelFva.global;
  if (!g || !g.length) { host.innerHTML = '<p class="hint">FVA data unavailable.</p>'; return; }
  const rxnName = (id) => (STATE.reactionsPanel && STATE.reactionsPanel[id] && STATE.reactionsPanel[id].name) || '';
  const top = g.slice(0, 20);
  const nmax = Math.max(...top.map((x) => x.n_forced_models), 1);
  const rows = top.map((x) =>
    `<div class="an-distrow" title="${escapeHtml(x.rxn)}${rxnName(x.rxn) ? ' — ' + rxnName(x.rxn) : ''}: flux-forced in ${x.n_forced_models} models, blocked in ${x.n_blocked_models}">` +
    `<span class="an-distlbl"><a href="#" class="rxn-link" data-rxn="${escapeHtml(x.rxn)}">${escapeHtml(x.rxn)}</a></span>` +
    `<span class="an-distbar"><span class="an-seg" style="width:${(x.n_forced_models / nmax * 100).toFixed(1)}%;background:var(--warn)"></span></span>` +
    `<span class="an-distn">${x.n_forced_models}</span></div>`).join('');
  host.innerHTML = `<div class="an-distwrap an-crit">${rows}</div>
    <p class="hint"># panel models in which each reaction is <strong>flux-forced</strong> (must carry flux for optimal growth) — obligate reactions of the growth backbone.</p>`;
  bindRxnLinks(host);
}

// -------------------- bootstrap --------------------
(async function init() {
  enhanceResizableGrids();
  await fetchHealth();
  await renderVariants();
  if (!STATE.staticMode) renderOverrides();
})();
