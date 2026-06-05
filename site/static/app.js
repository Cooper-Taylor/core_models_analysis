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
  selectedVariant: null,
  selectedRxn: null,
  rxnFilter: 'panel',
};

// -------------------- tab switching --------------------
document.querySelectorAll('nav button').forEach((btn) =>
  btn.addEventListener('click', () => {
    document.querySelectorAll('nav button').forEach((b) => b.classList.remove('active'));
    document.querySelectorAll('.tab').forEach((t) => t.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(`tab-${btn.dataset.tab}`).classList.add('active');
    if (btn.dataset.tab === 'reaction' && !STATE.reactionsPanel) loadReactions();
    if (btn.dataset.tab === 'sandbox') initSandbox();
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

async function renderVariants() {
  const m = await loadManifest();
  document.getElementById('manifest-meta').innerHTML =
    ` &nbsp;|&nbsp; ${m.variants.length} variants &nbsp;|&nbsp; ` +
    `${m.n_msdb_rxns.toLocaleString()} MSDB reactions &nbsp;|&nbsp; ` +
    `${m.n_panel_rxns.toLocaleString()} reactions in panel union or changed by ≥1 variant &nbsp;|&nbsp; ` +
    `built ${m.generated_at}`;

  const tbody = document.querySelector('#variants-table tbody');
  tbody.innerHTML = '';
  m.variants.forEach((v) => {
    const tr = document.createElement('tr');
    tr.dataset.tag = v.tag;
    tr.innerHTML = `
      <td><span class="tag">${escapeHtml(v.tag)}</span></td>
      <td>${escapeHtml(v.section)}</td>
      <td>${escapeHtml(v.title)}</td>
      <td class="num">${v.n_changed_vs_baseline.toLocaleString()}</td>
      <td class="num">${v.n_models_flux_change}</td>
      <td class="num">${v.n_models_flip}</td>
    `;
    tr.addEventListener('click', () => selectVariant(v.tag, tr));
    tbody.appendChild(tr);
  });
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
  renderVariantDetail(p);
}

function renderVariantDetail(p) {
  const pane = document.getElementById('variant-detail');
  let html = `
    <h3>${escapeHtml(p.tag)} — ${escapeHtml(p.title)}</h3>
    <p class="hint">${escapeHtml(p.section)}</p>
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
      </div>

      <h3>Panel FBA impact</h3>
      <div class="flux-impact-grid">
        <div class="card"><h4>panel size</h4><div class="stat">${p.panel_fba.length}</div></div>
        <div class="card"><h4>models grow-status flipped</h4><div class="stat">${p.n_models_flip}</div></div>
        <div class="card"><h4>models with flux Δ &gt; 1e-6</h4><div class="stat">${p.n_models_flux_change}</div></div>
      </div>

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
      ${p.diffs.length > 50 ? `<p class="hint">… ${p.diffs.length - 50} more not shown. Use the Reaction Explorer to browse.</p>` : ''}

      <h3>Models that changed flux (top 25 by |Δ|)</h3>
      <table class="changed-by-table">
        <thead><tr><th>model_id</th><th>baseline grows</th><th>variant grows</th><th>baseline flux</th><th>variant flux</th><th>Δ flux</th></tr></thead>
        <tbody>
          ${[...p.panel_fba].filter((r) => Math.abs(r.delta_flux) > 1e-6)
            .sort((a, b) => Math.abs(b.delta_flux) - Math.abs(a.delta_flux))
            .slice(0, 25).map((r) =>
              `<tr><td>${escapeHtml(r.model_id)}</td>
                   <td>${r.baseline_grows ? '✓' : '✗'}</td>
                   <td>${r.variant_grows ? '✓' : '✗'}</td>
                   <td class="num">${r.baseline_flux.toFixed(4)}</td>
                   <td class="num">${r.variant_flux.toFixed(4)}</td>
                   <td class="num ${r.delta_flux > 0 ? 'diff-up' : 'diff-down'}">${(r.delta_flux >= 0 ? '+' : '') + r.delta_flux.toFixed(4)}</td></tr>`
            ).join('')}
        </tbody>
      </table>`;
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
  renderReactionList();
}

async function ensureReactionsOther() {
  if (!STATE.reactionsOther) {
    STATE.reactionsOther = await API.data('reactions_other.json');
  }
  return STATE.reactionsOther;
}

function renderReactionList() {
  const q = document.getElementById('rxn-search').value.toLowerCase().trim();
  const filt = document.getElementById('rxn-filter').value;
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
  if (q) {
    entries = entries.filter((r) =>
      (r.id || '').toLowerCase().includes(q) ||
      (r.name || '').toLowerCase().includes(q) ||
      (r.definition || '').toLowerCase().includes(q));
  }
  entries.sort((a, b) => (a.id || '').localeCompare(b.id || ''));
  const ul = document.getElementById('reaction-list');
  document.getElementById('rxn-result-count').textContent = `${entries.length.toLocaleString()} reactions`;
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
  if (!r) {
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

    ${stoich ? `<h3>stoichiometry</h3><ul class="stoich-list">${stoich}</ul>` : ''}

    ${changedRows ? `<h3>Changed by ${r.changed_by.length} variant${r.changed_by.length > 1 ? 's' : ''}</h3>
      <table class="changed-by-table">
        <thead><tr><th>variant</th><th>baseline</th><th>variant</th></tr></thead>
        <tbody>${changedRows}</tbody>
      </table>` : `<p class="hint">No variants change this reaction's direction vs baseline.</p>`}

    <h3>Per-mode panel sweep <span class="hint">— live FBA</span></h3>
    <p class="hint">Run FBA on panel models that contain this reaction under three forced modes:
      <strong>off</strong> (lb=ub=0), <strong>forward</strong> (lb=0, ub=1000), <strong>reverse</strong> (lb=−1000, ub=0).
      Comparison is vs the unaltered cascade (baseline by default).
      First call to a model loads the JSON; expect ~5-20s.</p>
    <div class="run-row">
      <label>variant:
        <select id="rxn-sweep-variant"></select>
      </label>
      <button id="rxn-sweep-run" class="primary">Run sweep</button>
      <span id="rxn-sweep-status" class="hint"></span>
    </div>
    <div id="rxn-sweep-results"></div>
  `;

  // Populate variant select.
  loadManifest().then((m) => {
    const sel = document.getElementById('rxn-sweep-variant');
    sel.innerHTML = '';
    m.variants.forEach((v) => {
      const o = document.createElement('option');
      o.value = v.tag;
      o.textContent = `${v.tag} — ${v.title}`;
      sel.appendChild(o);
    });
  });

  document.getElementById('rxn-sweep-run').addEventListener('click', async () => {
    const variant = document.getElementById('rxn-sweep-variant').value;
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
        variant,
        modes: ['as_is', 'off', 'forward', 'reverse', 'free'],
      });
      status.textContent = `done in ${(performance.now() - t0) / 1000 | 0}s ` +
        `(${res.n_models} models, ${(res.elapsed_s).toFixed(1)}s server-side)`;
      renderSweep(res, out);
    } catch (exc) {
      status.textContent = 'error: ' + exc.message;
    } finally {
      btn.disabled = false;
    }
  });
}

function renderSweep(res, out) {
  const modes = Object.keys(res.by_mode);
  const cards = modes.map((mode) => {
    const rows = Object.entries(res.by_mode[mode]);
    const grew = rows.filter(([_, v]) => v.grows).length;
    const flipped = rows.filter(([mid, v]) => v.grows !== (res.baseline[mid] || {}).grows).length;
    const fluxChanged = rows.filter(([mid, v]) => Math.abs(v.delta_flux) > 1e-6).length;
    return `<div class="card">
      <h4>mode: ${mode}</h4>
      <div class="stat">${grew}/${rows.length}<span class="stat-sub"> grow</span></div>
      <div class="stat" style="font-size:14px">${flipped} flipped, ${fluxChanged} flux Δ</div>
    </div>`;
  }).join('');

  // Per-model table
  const mids = Object.keys(res.baseline).sort();
  const rows = mids.map((mid) => {
    const baseG = res.baseline[mid].grows;
    const baseF = res.baseline[mid].growth_flux;
    const cells = modes.map((mode) => {
      const v = res.by_mode[mode][mid];
      if (!v) return '<td>—</td>';
      const flipped = v.grows !== baseG;
      const flag = flipped ? (v.grows ? ' ↑' : ' ↓') : '';
      const cls = v.delta_flux > 0 ? 'diff-up' : (v.delta_flux < 0 ? 'diff-down' : 'diff-zero');
      return `<td class="num ${cls}">${v.growth_flux.toFixed(4)}${flag}</td>`;
    }).join('');
    return `<tr><td>${escapeHtml(mid)}</td><td class="num">${baseF.toFixed(4)}</td>${cells}</tr>`;
  }).join('');

  out.innerHTML = `
    <div class="flux-impact-grid">${cards}</div>
    <h4>Per-model growth flux</h4>
    <table class="changed-by-table">
      <thead><tr><th>model_id</th><th>baseline flux</th>${modes.map((m) => `<th>${m}</th>`).join('')}</tr></thead>
      <tbody>${rows}</tbody>
    </table>
    <p class="hint">"↑" = model became a grower under the mode; "↓" = model stopped growing. Color = sign of Δ flux vs the chosen variant's baseline.</p>
  `;
}

// Search wiring
document.getElementById('rxn-search').addEventListener('input', () => renderReactionList());
document.getElementById('rxn-filter').addEventListener('change', (e) => {
  STATE.rxnFilter = e.target.value;
  if (STATE.rxnFilter !== 'panel') ensureReactionsOther().then(renderReactionList);
  else renderReactionList();
});

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
  const rxn = document.getElementById('ov-rxn').value.trim();
  const mode = document.getElementById('ov-mode').value;
  if (!rxn) return;
  SANDBOX_OVERRIDES.push({ rxn, mode });
  document.getElementById('ov-rxn').value = '';
  renderOverrides();
});

document.getElementById('sandbox-run').addEventListener('click', async () => {
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
  await renderVariants();
  renderOverrides();
})();
