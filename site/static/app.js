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
  selectedVariant: null,
  selectedRxn: null,
  rxnFilter: 'panel',
  rxnVariantFilter: 'any',     // 'any' / 'none' / specific tag
  rxnFluxImpactedOnly: false,  // subfilter for §3c
  scope: 'panel',              // 'panel' | 'all'
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

// -------------------- tab switching --------------------
document.querySelectorAll('nav button').forEach((btn) =>
  btn.addEventListener('click', () => {
    document.querySelectorAll('nav button').forEach((b) => b.classList.remove('active'));
    document.querySelectorAll('.tab').forEach((t) => t.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(`tab-${btn.dataset.tab}`).classList.add('active');
    if (btn.dataset.tab === 'reaction' && !STATE.reactionsPanel) loadReactions();
    if (btn.dataset.tab === 'sandbox' && !STATE.staticMode) initSandbox();
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
      ? `${fmt(inst.count)} <span class="hint" style="font-size:11px">(${fmt(inst.n_models_with_change)} models)</span>`
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
      `baseline: ${r.baseline_flux.toFixed(4)} (${r.baseline_grows ? 'grew' : 'no growth'})\n` +
      `variant:  ${r.variant_flux.toFixed(4)} (${r.variant_grows ? 'grew' : 'no growth'})\n` +
      `Δ flux:   ${r.delta_flux.toFixed(4)}${flipped ? '   ← grow-status flipped' : ''}`;
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
        bar scale: ±${max_abs.toFixed(2)} flux units
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
        <div class="stat-sub"><span class="diff-up">${summary.n_not_default_now_grew.toLocaleString()}</span> didn't grow → now do</div>
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
      <thead><tr><th>model_id</th><th>baseline flux</th><th>variant flux</th><th>Δ flux</th><th>grow base → variant</th></tr></thead>
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
      <div class="transition-grid">
        ${Object.entries(p.transitions || {}).sort().map(([k, v]) => {
          const [from, to] = k.split('->');
          const dir = (to === '>' || to === '<') && (from === '=' || from === '?') ? 'up'
                      : (from === '>' || from === '<') && (to === '=' || to === '?') ? 'down' : '';
          return `<div class="t-cell ${dir}">${revBadge(from)} → ${revBadge(to)}<span>${v}</span></div>`;
        }).join('')}
      </div>`;

    if (STATE.scope === 'panel') {
      html += `
        <h3>Panel FBA impact</h3>
        <div class="flux-impact-grid">
          <div class="card"><h4>panel size</h4><div class="stat">${p.panel_fba.length}</div></div>
          <div class="card"><h4>models grow-status flipped</h4><div class="stat">${p.n_models_flip}</div></div>
          <div class="card"><h4>models with flux Δ &gt; 1e-6</h4><div class="stat">${p.n_models_flux_change}</div></div>
        </div>`;
    } else {
      const sm = STATE.allModelsSummary?.variants?.[p.tag];
      html += `
        <h3>All-models FBA impact <span class="hint">— ${STATE.allModelsSummary?.n_all_models?.toLocaleString() || '?'} core models</span></h3>
        ${renderAllModelsStatsCard(p, allFbaRows)}`;
    }

    html += `
      <h3>Top reactions changed (first 50)</h3>
      <table class="changed-by-table">
        <thead><tr><th>rxn</th><th>base</th><th>new</th><th></th></tr></thead>
        <tbody>
          ${p.diffs.slice(0, 50).map((d) =>
            `<tr><td><a href="#" class="rxn-link" data-rxn="${escapeHtml(d.rxn)}">${escapeHtml(d.rxn)}</a></td>
                 <td>${revBadge(d.base)}</td>
                 <td>${revBadge(d.new)}</td>
                 <td></td></tr>`
          ).join('')}
        </tbody>
      </table>
      ${p.diffs.length > 50 ? `<p class="hint">… ${p.diffs.length - 50} more not shown. Use the Reaction Explorer to browse.</p>` : ''}`;

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

  // Cross-link rxn IDs in the variant view → reaction explorer
  pane.querySelectorAll('.rxn-link').forEach((a) =>
    a.addEventListener('click', (e) => {
      e.preventDefault();
      document.querySelector('nav button[data-tab="reaction"]').click();
      setTimeout(() => selectRxn(a.dataset.rxn), 30);
    }));
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
  renderReactionList();
}

// Cache: tag -> Set of rxn_ids whose direction-change has a downstream
// effect on biomass flux in at least one panel model containing it.
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
      ${r.pathways && r.pathways.length ? `<dt>pathways</dt><dd>${r.pathways.map(escapeHtml).join(', ')}</dd>` : ''}
      ${r.in_panel !== undefined ? `<dt>in panel</dt><dd>${r.in_panel ? `yes (${r.panel_freq} of 100 models)` : 'no'}</dd>` : ''}
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
  const rowsHtml = mids.map((mid) => {
    const b = res.baseline[mid];
    const v = modeRows[mid];
    const baseF = b.growth_flux;
    if (b.grows) nGrewBefore += 1;
    if (!v) {
      return `<tr><td><code>${escapeHtml(mid)}</code></td>
        <td class="num">${baseF.toFixed(4)}${b.grows ? '' : ' ✗'}</td>
        <td class="num">—</td>
        <td class="num">—</td></tr>`;
    }
    const newF = v.growth_flux;
    if (v.grows) nGrewAfter += 1;
    const flipped = v.grows !== b.grows;
    if (flipped) nFlipped += 1;
    if (Math.abs(v.delta_flux) > 1e-6) nFluxChanged += 1;
    const flag = flipped ? (v.grows ? ' ↑' : ' ↓') : '';
    const cls = v.delta_flux > 0 ? 'diff-up' : (v.delta_flux < 0 ? 'diff-down' : 'diff-zero');
    return `<tr${flipped ? ' class="flipped"' : ''}>
      <td><code>${escapeHtml(mid)}</code></td>
      <td class="num">${baseF.toFixed(4)}${b.grows ? '' : ' ✗'}</td>
      <td class="num ${cls}">${newF.toFixed(4)}${flag}</td>
      <td class="num ${cls}">${v.delta_flux >= 0 ? '+' : ''}${v.delta_flux.toFixed(4)}</td>
    </tr>`;
  }).join('');

  out.innerHTML = `
    <div class="flux-impact-grid">
      <div class="card">
        <h4>variant</h4>
        <div class="stat" style="font-size:14px"><code>${escapeHtml(variant)}</code></div>
        <div class="stat-sub">new direction = ${escapeHtml(newRev)} (mode: ${escapeHtml(mode)})</div>
      </div>
      <div class="card">
        <h4>panel grow-status</h4>
        <div class="stat">${nGrewBefore} → ${nGrewAfter}</div>
        <div class="stat-sub">${nFlipped} flipped, ${nFluxChanged} flux Δ</div>
      </div>
      <div class="card">
        <h4>panel models</h4>
        <div class="stat">${mids.length}</div>
        <div class="stat-sub">containing this reaction</div>
      </div>
    </div>
    <h4>Per-model growth flux</h4>
    <table class="changed-by-table">
      <thead><tr>
        <th>model_id</th>
        <th>default flux</th>
        <th>variant flux (${escapeHtml(newRev)})</th>
        <th>Δ flux</th>
      </tr></thead>
      <tbody>${rowsHtml}</tbody>
    </table>
    <p class="hint">"↑" = model became a grower under the variant's new direction for this reaction; "↓" = model stopped growing. "✗" = was not a grower in the baseline.</p>
  `;
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
  // Disable the flux-impacted-only subfilter when "any/none" is selected
  // (it only makes sense for a specific variant).
  const sub = document.getElementById('rxn-flux-impacted-only');
  const subEnabled = !(e.target.value === 'any' || e.target.value === 'none');
  sub.disabled = !subEnabled;
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
        <div class="badge">mean flux: <span class="n">${mean.toFixed(4)}</span></div>
      </div>
    </div>
    <table class="changed-by-table">
      <thead><tr><th>model_id</th><th>status</th><th>grows</th><th>growth_flux</th><th>n overrides applied</th></tr></thead>
      <tbody>${rs.sort((a, b) => b.growth_flux - a.growth_flux).map((r) =>
        `<tr><td>${escapeHtml(r.model_id)}</td>
              <td>${escapeHtml(r.status)}</td>
              <td>${r.grows ? '✓' : '✗'}</td>
              <td class="num">${r.growth_flux.toFixed(4)}</td>
              <td class="num">${r.n_overrides}</td></tr>`
      ).join('')}</tbody>
    </table>`;
}

// -------------------- bootstrap --------------------
(async function init() {
  await fetchHealth();
  await renderVariants();
  if (!STATE.staticMode) renderOverrides();
})();
