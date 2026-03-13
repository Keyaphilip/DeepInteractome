/* ── DeepInteractome Web UI — app.js ──────────────────────────────────────── */

const API_BASE = '';   // same origin; update if hosting separately

// ── Particle system ───────────────────────────────────────────────────────────
(function spawnParticles() {
  const container = document.getElementById('heroParticles');
  if (!container) return;
  for (let i = 0; i < 28; i++) {
    const p = document.createElement('div');
    p.className = 'particle';
    const size = Math.random() * 4 + 2;
    p.style.cssText = `
      width: ${size}px; height: ${size}px;
      left: ${Math.random() * 100}%;
      bottom: ${Math.random() * -20}%;
      animation-duration: ${Math.random() * 14 + 10}s;
      animation-delay: ${Math.random() * -15}s;
    `;
    container.appendChild(p);
  }
})();

// ── Variant row template ───────────────────────────────────────────────────────
let rowCount = 0;

function createVariantRow(defaults = {}) {
  rowCount++;
  const id = rowCount;
  const div = document.createElement('div');
  div.className = 'variant-row';
  div.dataset.rowId = id;
  div.innerHTML = `
    <div class="field-group">
      <label class="field-label">CHROM</label>
      <input class="field-input" id="chrom_${id}" placeholder="chr17"
             value="${defaults.chrom || ''}" />
    </div>
    <div class="field-group">
      <label class="field-label">POS</label>
      <input class="field-input" id="pos_${id}" type="number" placeholder="7673802"
             value="${defaults.pos || ''}" />
    </div>
    <div class="field-group">
      <label class="field-label">REF</label>
      <input class="field-input" id="ref_${id}" placeholder="C" maxlength="50"
             value="${defaults.ref || ''}" />
    </div>
    <div class="field-group">
      <label class="field-label">ALT</label>
      <input class="field-input" id="alt_${id}" placeholder="T" maxlength="50"
             value="${defaults.alt || ''}" />
    </div>
    <div class="field-group">
      <label class="field-label">AF (0–1)</label>
      <input class="field-input" id="af_${id}" type="number" step="0.0001"
             min="0" max="1" placeholder="0.001"
             value="${defaults.af !== undefined ? defaults.af : ''}" />
    </div>
    <button class="remove-btn" title="Remove variant" data-row="${id}">✕</button>
  `;
  div.querySelector('.remove-btn').addEventListener('click', () => removeRow(id));
  return div;
}

function addRow(defaults = {}) {
  const container = document.getElementById('variantRows');
  container.appendChild(createVariantRow(defaults));
  updateRemoveButtons();
}

function removeRow(id) {
  const row = document.querySelector(`[data-row-id="${id}"]`);
  if (row) row.remove();
  updateRemoveButtons();
}

function updateRemoveButtons() {
  const rows = document.querySelectorAll('.variant-row');
  rows.forEach(r => {
    const btn = r.querySelector('.remove-btn');
    btn.style.visibility = rows.length > 1 ? 'visible' : 'hidden';
  });
}

function collectVariants() {
  const rows = document.querySelectorAll('.variant-row');
  const variants = [];
  const errors = [];

  rows.forEach((row, idx) => {
    const id = row.dataset.rowId;
    const chrom = document.getElementById(`chrom_${id}`).value.trim();
    const pos   = parseInt(document.getElementById(`pos_${id}`).value, 10);
    const ref   = document.getElementById(`ref_${id}`).value.trim().toUpperCase();
    const alt   = document.getElementById(`alt_${id}`).value.trim().toUpperCase();
    const afRaw = document.getElementById(`af_${id}`).value.trim();
    const af    = afRaw === '' ? 0.0 : parseFloat(afRaw);

    if (!chrom) errors.push(`Row ${idx + 1}: CHROM is required`);
    if (!pos || isNaN(pos)) errors.push(`Row ${idx + 1}: POS must be a number`);
    if (!ref) errors.push(`Row ${idx + 1}: REF is required`);
    if (!alt) errors.push(`Row ${idx + 1}: ALT is required`);
    if (isNaN(af) || af < 0 || af > 1) errors.push(`Row ${idx + 1}: AF must be between 0 and 1`);

    if (!errors.length) variants.push({ chrom, pos, ref, alt, af });
  });

  return { variants, errors };
}

// ── Example data ──────────────────────────────────────────────────────────────
const EXAMPLES = [
  { chrom: 'chr17', pos: 7673802,  ref: 'C',   alt: 'T',  af: 0.0001  },
  { chrom: 'chr13', pos: 32930644, ref: 'A',   alt: 'T',  af: 0.0003  },
  { chrom: 'chr1',  pos: 10177,    ref: 'A',   alt: 'AC', af: 0.425   },
];

document.getElementById('loadExampleBtn').addEventListener('click', () => {
  document.getElementById('variantRows').innerHTML = '';
  rowCount = 0;
  EXAMPLES.forEach(addRow);
});

// ── CSV Upload ─────────────────────────────────────────────────────────────────
document.getElementById('csvUpload').addEventListener('change', (e) => {
  const file = e.target.files[0];
  if (!file) return;

  const reader = new FileReader();
  reader.onload = (event) => {
    const text = event.target.result;
    const lines = text.split(/\r?\n/).filter(l => l.trim().length > 0);
    if (lines.length < 2) {
      alert("CSV must contain a header row and at least one data row.");
      return;
    }
    
    // Naive CSV parsing
    const headers = lines[0].split(',').map(h => h.trim().toLowerCase());
    const chromIdx = headers.indexOf('chrom');
    const posIdx = headers.indexOf('pos');
    const refIdx = headers.indexOf('ref');
    const altIdx = headers.indexOf('alt');
    const afIdx = headers.indexOf('af');

    if (chromIdx === -1 || posIdx === -1 || refIdx === -1 || altIdx === -1) {
      alert("CSV must contain 'chrom', 'pos', 'ref', 'alt' columns.");
      return;
    }

    document.getElementById('variantRows').innerHTML = '';
    rowCount = 0;

    for (let i = 1; i < lines.length; i++) {
        // Limit to 50 rows to keep the UI smooth
        if (i > 50) {
            alert("Only loading the first 50 variants from the CSV to prevent UI lag.");
            break;
        }
        
        const row = lines[i].split(',').map(c => c.trim());
        const afVal = row[afIdx];
        const variant = {
            chrom: row[chromIdx],
            pos: parseInt(row[posIdx], 10),
            ref: row[refIdx],
            alt: row[altIdx],
            af: (afIdx !== -1 && afVal) ? parseFloat(afVal) : 0.0
        };
        addRow(variant);
    }
    
    e.target.value = ""; // Clear file input
  };
  reader.readAsText(file);
});

// ── Form submit ────────────────────────────────────────────────────────────────
document.getElementById('variantForm').addEventListener('submit', async (e) => {
  e.preventDefault();
  const { variants, errors } = collectVariants();

  const errEl   = document.getElementById('resultsError');
  const tableWrap = document.getElementById('resultsTableWrap');
  const placeholder = document.getElementById('resultsPlaceholder');
  const statsBar  = document.getElementById('statsBar');

  // Clear previous
  errEl.classList.add('hidden');
  tableWrap.classList.add('hidden');
  placeholder.classList.remove('hidden');
  statsBar.style.display = 'none';

  if (errors.length) {
    errEl.textContent = errors.join(' · ');
    errEl.classList.remove('hidden');
    return;
  }

  // Loading state
  const btn = document.getElementById('predictBtn');
  const btnText = btn.querySelector('.btn-text');
  const spinner = document.getElementById('spinner');
  btn.disabled = true;
  btnText.textContent = 'Predicting…';
  spinner.classList.remove('hidden');

  try {
    const resp = await fetch(`${API_BASE}/predict`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ variants }),
    });

    if (!resp.ok) {
      const text = await resp.text();
      throw new Error(`API error ${resp.status}: ${text}`);
    }

    const data = await resp.json();
    renderResults(data);

  } catch (err) {
    errEl.textContent = `Request failed: ${err.message}`;
    errEl.classList.remove('hidden');
    placeholder.classList.remove('hidden');
  } finally {
    btn.disabled = false;
    btnText.textContent = 'Predict Pathogenicity';
    spinner.classList.add('hidden');
  }
});

// ── Render results ─────────────────────────────────────────────────────────────
function renderResults(data) {
  const tbody = document.getElementById('resultsBody');
  const placeholder = document.getElementById('resultsPlaceholder');
  const tableWrap = document.getElementById('resultsTableWrap');
  const modelBadge = document.getElementById('modelBadge');
  const statsBar = document.getElementById('statsBar');

  tbody.innerHTML = '';
  placeholder.classList.add('hidden');
  tableWrap.classList.remove('hidden');

  // Model badge
  modelBadge.textContent = data.model_used || 'unknown';
  modelBadge.classList.remove('hidden');

  let pathCount = 0, benignCount = 0, totalProb = 0;

  data.predictions.forEach((pred, i) => {
    const isPath = pred.result === 'PATHOGENIC';
    if (isPath) pathCount++; else benignCount++;
    totalProb += pred.pathogenic_probability;

    const pct = (pred.pathogenic_probability * 100).toFixed(1);
    const fillClass = isPath ? 'fill-pathogenic' : 'fill-benign';
    const badgeClass = isPath ? 'badge-pathogenic' : 'badge-benign';
    const emoji = isPath ? '🔴' : '🟢';

    const tr = document.createElement('tr');
    tr.style.animationDelay = `${i * 60}ms`;
    tr.innerHTML = `
      <td class="variant-id">
        ${esc(pred.chrom)}:${pred.pos} ${esc(pred.ref)}→${esc(pred.alt)}
      </td>
      <td>
        <span class="badge ${badgeClass}">${emoji} ${pred.result}</span>
      </td>
      <td>
        <div class="confidence-bar-wrap">
          <div class="confidence-bar-bg">
            <div class="confidence-bar-fill ${fillClass}"
                 style="width: ${pct}%"></div>
          </div>
          <span class="confidence-pct">${pct}%</span>
        </div>
      </td>
    `;
    tbody.appendChild(tr);
  });

  // Stats
  const total = data.predictions.length;
  document.getElementById('statTotal').textContent  = total;
  document.getElementById('statPath').textContent   = pathCount;
  document.getElementById('statBenign').textContent = benignCount;
  document.getElementById('statAvg').textContent    = `${(totalProb / total * 100).toFixed(1)}%`;
  statsBar.style.display = 'flex';
}

function esc(str) {
  return String(str)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;');
}

// ── Init ───────────────────────────────────────────────────────────────────────
document.getElementById('addRowBtn').addEventListener('click', () => addRow());
addRow();  // default first row
