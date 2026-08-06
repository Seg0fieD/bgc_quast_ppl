/* global Chart */

/* ---------------------------------------------------------------------------
 * BGC-QUAST RUN SUMMARY
 * ------------------------------------------------------------------------- */

function renderRunSummary(meta) {
    const el = document.getElementById('runSummary');
    if (!el) return;

    const requested = meta.requested_mode;
    const resolved  = meta.running_mode;

    let modeText;
    if (requested === 'auto') {
        modeText = `The running mode is <code>auto</code> that determined the mode as <code>${resolved}</code>.`;
    } else {
        modeText = `The running mode is <code>${resolved}</code>.`;
    }

    let text = `${modeText} All statistics are based on BGCs of size ≥ <code>${meta.min_bgc_length}</code> bp. `;

    const margin = meta.bgc_completeness_margin;
    if (margin != null) {
        text += `BGCs are considered complete if located at least <code>${margin}</code> bp from contig edges.`;
    }

    if (resolved === 'compare_tools' && meta.compare_tools_overlap_threshold != null) {
        const pct = Math.round(meta.compare_tools_overlap_threshold * 100);
        text += ` The BGC overlap threshold is set to <code>${pct}%</code>.`;
    }

    el.innerHTML = text;
}


/* ---------------------------------------------------------------------------
 * HEATMAP UTILITIES
 * ------------------------------------------------------------------------- */

/**
 * Compute the median of a sorted numeric array
 * @param {number[]} arr
 *  @returns {number}
 */
function getMedian(arr) {
    const mid = Math.floor(arr.length / 2); // Calculate the middle index of the array
    return arr.length % 2 === 0
        ? (arr[mid - 1] + arr[mid]) / 2
        : arr[mid];
}

// Color helpers: hex <-> rgb + mixing
// Convert hex → RGB
function hexToRgb(hex) {
    const c = hex.replace('#', '');
    return {
        r: parseInt(c.substring(0, 2), 16),
        g: parseInt(c.substring(2, 4), 16),
        b: parseInt(c.substring(4, 6), 16)
    };
}

// Convert RGB → hex
function rgbToHex({ r, g, b }) {
    const toHex = x => x.toString(16).padStart(2, '0');
    return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}

/**
 * Mix a base color with white.
 * factor = 0   → white
 * factor = 1   → base color
 */
function mixColorWithWhite(baseHex, factor) {
    const base = hexToRgb(baseHex);
    const white = { r: 255, g: 255, b: 255 };

    return rgbToHex({
        r: Math.round(white.r + (base.r - white.r) * factor),
        g: Math.round(white.g + (base.g - white.g) * factor),
        b: Math.round(white.b + (base.b - white.b) * factor)
    });
}

/**
 * Wrapper: returns a hex color that is a lighter/darker
 * version of baseHex, controlled by factor [0..1].
 */
function getColorHex(baseHex, factor) {
    return mixColorWithWhite(baseHex, factor);
}

// Base colors for the heatmap:
const LOW_COLOR_HEX  = "#80852a";  // olive green-yellow
const HIGH_COLOR_HEX = "#396B9E";  // denim blue

const FACT_WHITE = 0.0;   // ~100% lightness (white)
const FACT_MIN   = 0.3;   // ~75% lightness
const FACT_INNER = 0.55;  // ~65% lightness
const FACT_OUTER = 0.8;   // ~55% lightness

/**
 * Draw the heatmap legend (smallest → median → largest) above the table.
 * Uses LOW_COLOR_HEX and HIGH_COLOR_HEX.
 */
function drawHeatmapLegend() {
    const canvas = document.getElementById('gradientHeatmap');
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    const gradient = ctx.createLinearGradient(0, 0, width, 0);

    // Smallest → median → largest
    gradient.addColorStop(0,   getColorHex(LOW_COLOR_HEX, FACT_INNER));
    gradient.addColorStop(0.5, "#FFFFFF"); // median
    gradient.addColorStop(1,   getColorHex(HIGH_COLOR_HEX, FACT_INNER));

    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, width, height);
}

/**
 * Apply background color to cells in a single row based on outliers/statistics.
 * @param {HTMLTableCellElement[]} cells
 * @param {number[]} values
 * @param {'more_is_better'|'less_is_better'} [direction='more_is_better']
 */
function heatMapOneRow(cells, values, direction = 'more_is_better') {
    // Skip if too few values or all are the same (no variation to color)
    if (values.length < 2 || new Set(values).size === 1) return;

    // Sort values and calculate statistical parameters
    const sorted = [...values].sort((a, b) => a - b);
    const median = getMedian(sorted);
    const q1 = sorted[Math.floor((sorted.length - 1) / 4)];
    const q3 = sorted[Math.floor((sorted.length - 1) * 3 / 4)];
    const iqr = q3 - q1;

    // Calculate "fences" for mild/extreme outliers
    const lowOuter = q1 - 3 * iqr;
    const lowInner = q1 - 1.5 * iqr;
    const topInner = q3 + 1.5 * iqr;
    const topOuter = q3 + 3 * iqr;


    // Choose which side is "low" and which is "high" based on metric direction
    const [lowColorHex, topColorHex] =
        direction === 'less_is_better'
            ? [HIGH_COLOR_HEX, LOW_COLOR_HEX]  // reversed
            : [LOW_COLOR_HEX,  HIGH_COLOR_HEX];

    // Apply coloring logic per cell
    for (let i = 0; i < cells.length; i++) {
        const cell = cells[i];
        const num = values[i];
        if (num == null || Number.isNaN(num)) continue;

        if (num < lowOuter) {
            // Extreme low outlier
            cell.style.backgroundColor = getColorHex(lowColorHex, FACT_OUTER);
            cell.style.color = 'white';
        } else if (num < lowInner) {
            // Mild low outlier
            cell.style.backgroundColor = getColorHex(lowColorHex, FACT_INNER);
        } else if (num < median) {
            // Between lowInner and median → interpolate white → lowColor
            const denom = (median - lowInner) || 1;
            const k = (FACT_MIN - FACT_WHITE) / denom;
            const factor = FACT_WHITE + (median - num) * k;
            cell.style.backgroundColor = getColorHex(lowColorHex, factor);
        } else if (num > topOuter) {
            // Extreme high outlier
            cell.style.backgroundColor = getColorHex(topColorHex, FACT_OUTER);
            cell.style.color = 'white';
        } else if (num > topInner) {
            // Mild high outlier
            cell.style.backgroundColor = getColorHex(topColorHex, FACT_INNER);
        } else if (num > median) {
            // Between median and topInner → interpolate white → topColor
            const denom = (topInner - median) || 1;
            const k = (FACT_MIN - FACT_WHITE) / denom;
            const factor = FACT_WHITE + (num - median) * k;
            cell.style.backgroundColor = getColorHex(topColorHex, factor);
        }
        // num === median → stays white
    }
}

/**
 * Global store to re-apply/remove heatmap colors when toggling the heatmap.
 * Each entry: { cells: HTMLTableCellElement[], values: number[] }
 */
const allNumericCells = [];
let showProductRows = false;
let showCompletenessRows = false;

function classifyRowLabel(rawLabel) {
    const label = String(rawLabel ?? '').trim();
    const m = /^(.*?)\s*\((.*)\)$/.exec(label);

    if (!m) {
        return { isTotal: true, isProductOnly: false, isCompletenessOnly: false, isFull: false };
    }

    const insideParts = m[2].split(',').map(s => s.trim().toLowerCase());
    const hasCompleteness = insideParts.some(v => ['complete', 'incomplete', 'unknown completeness'].includes(v));

    if (insideParts.length === 1) {
        if (hasCompleteness) {
            return { isTotal: false, isProductOnly: false, isCompletenessOnly: true, isFull: false };
        }
        return { isTotal: false, isProductOnly: true, isCompletenessOnly: false, isFull: false };
    }

    return { isTotal: false, isProductOnly: false, isCompletenessOnly: false, isFull: true };
}

function applyExtendedState() {
    const rows = document.querySelectorAll('#reportTableContainer table tr');

    rows.forEach((row, idx) => {
        if (idx === 0) return; // header row

        const headerCell = row.querySelector('th');
        if (!headerCell) return;

        const rawLabel = headerCell.textContent || '';
        if (rawLabel === 'Genome mining tool') {
            row.style.display = 'table-row';
            return;
        }

        const cls = classifyRowLabel(rawLabel);

        if (cls.isTotal) {
            row.style.display = 'table-row';
        } else if (cls.isProductOnly) {
            row.style.display = showProductRows ? 'table-row' : 'none';
        } else if (cls.isCompletenessOnly) {
            row.style.display = showCompletenessRows ? 'table-row' : 'none';
        } else if (cls.isFull) {
            row.style.display = (showProductRows && showCompletenessRows) ? 'table-row' : 'none';
        }
    });

    const productBtn = document.getElementById('toggleProductBtn');
    if (productBtn) {
        productBtn.textContent = showProductRows
            ? 'Hide stratification by product'
            : 'Show stratification by product';
        productBtn.classList.toggle('is-collapsed', showProductRows);
    }

    const completenessBtn = document.getElementById('toggleCompletenessBtn');
    if (completenessBtn) {
        completenessBtn.textContent = showCompletenessRows
            ? 'Hide stratification by completeness'
            : 'Show stratification by completeness';
        completenessBtn.classList.toggle('is-collapsed', showCompletenessRows);
    }

    const layout = document.getElementById('reportLayout');
    if (layout) {
        if (showProductRows || showCompletenessRows) layout.classList.add('extended-locked');
        else layout.classList.remove('extended-locked');
    }
}

/* ---------------------------------------------------------------------------
 * TABLE BUILDER (main report table + extended rows)
 * ------------------------------------------------------------------------- */
const isReferenceMode = (typeof reportMode !== 'undefined' && reportMode === "compare_to_reference");
let includeReferenceInHeatmap = false;
/**
 * Build the main report table and insert it into #reportTableContainer.
 * @param {Array[]} data - 2D array representing the report pivot table.
 */
function buildTable(data) {
    allNumericCells.length = 0;
    const container = document.getElementById('reportTableContainer');
    const table = document.createElement('table');
    const headerRow = document.createElement('tr');

    // Build column headers
    headerRow.appendChild(document.createElement('th')); // top-left empty cell
    data[0].slice(1).forEach((col, idx) => {
        const th = document.createElement('th');
        th.textContent = col;

        // idx=0 corresponds to j===1 (first data column) -> reference in reference mode
        if (isReferenceMode && idx === 0) {
            th.classList.add('ref-col-header');
        }

        headerRow.appendChild(th);
    });
    headerRow.classList.add('column-header-row');
    table.appendChild(headerRow);

    // Loop over each data row
    for (let i = 1; i < data.length; i++) {
        const row = document.createElement('tr');
        const rowData = data[i];

        // Hide extended rows by default (only show "total" by default)
        const label = String(rowData[0] ?? '');
        const isMiningTool = label.trim().toLowerCase() === 'genome mining tool';
        if (!isMiningTool) {
            row.classList.add('extended-row');
        }


        const numericCells = [];
        const numericValues = [];

        // Build each cell in the row
        data[i].forEach((cell, j) => {
            const isRowHeader = (j === 0);
            const td = document.createElement(isRowHeader ? 'th' : 'td');

            if (isRowHeader) {
                let labelText = (String(cell).toLowerCase() === 'genome mining tool')
                    ? 'Genome mining tool'
                    : String(cell ?? '');

                td.textContent = labelText;
                if (labelText === 'Genome mining tool') {
                    td.classList.add('row-label-total');
                }

                const rowClass = classifyRowLabel(labelText);

                if (rowClass.isTotal || labelText === 'Genome mining tool') {
                    td.classList.add('row-label-total');
                } else {
                    td.classList.add('row-label');
                }

                const rawLabel = String(cell ?? '');

                // Compare-tools mode-specific rows
                const isCompareToolsSpecific =
                    rawLabel.startsWith('# unique BGCs') ||
                    rawLabel.startsWith('Unique BGC rate');

                const isCompareRefSpecific =
                    rawLabel.startsWith('# full ref. BGCs') ||
                    rawLabel.startsWith('# full ref. BGCs, single-contig') ||
                    rawLabel.startsWith('# full ref. BGCs, multi-contig') ||
                    rawLabel.startsWith('# partial ref. BGCs') ||
                    rawLabel.startsWith('# missed ref. BGCs') ||
                    rawLabel.startsWith('# product type mismatches to ref.') ||
                    rawLabel.startsWith('# unmapped BGCs to ref.') ||
                    rawLabel.startsWith('Ref. BGC recovery rate');

                if ((reportMode === 'compare_tools' && isCompareToolsSpecific) ||
                    (reportMode === 'compare_to_reference' && isCompareRefSpecific)) {
                    td.classList.add('mode-specific-metric');
                }

            } else {
                td.textContent = cell;
            }

            // Only collect numeric values for heatmap
            if (j > 0) {
                const num = parseFloat(cell);

                // Only consider numeric cells for heatmap
                if (!Number.isNaN(num)) {
                    // In COMPARE_TO_REFERENCE mode, the first data column (j === 1)
                    // is the reference column
                    const isReferenceColumn = isReferenceMode && (j === 1);

                    if (isReferenceColumn && !includeReferenceInHeatmap) {
                        // Reference column stays uncolored (white) when heatmap is on.
                    } else {
                        // Normal numeric cell, or reference when included:
                        // participate in the heatmap.
                        numericCells.push(td);
                        numericValues.push(num);
                    }
                }
            }

            row.appendChild(td);
        });

        // Apply the statistical heatmap coloring
        heatMapOneRow(numericCells, numericValues);
        allNumericCells.push({ cells: numericCells, values: numericValues });

        table.appendChild(row);
    }

    container.appendChild(table);
}

/* ---------------------------------------------------------------------------
 * BAR CHART BUILDER (Chart.js) – totals / by type / by completeness
 * ------------------------------------------------------------------------- */

//product colors (modifiable
const productColors = {
    'NRPS': '#2e8b57',
    'PKS': '#f4a460',
    'RiPP': '#4169e1',
    'hybrid': '#689cba',
    'terpene': '#800080',
    'saccharide': '#f0c107',
    'alkaloid': '#dda0dd',
    'other': '#b56576',
    'unknown': '#000000',
};

/**
 * Lighten a HEX color toward white by a factor in [0, 1].
 * @param {string} hex
 * @param {number} [factor=0.5]
 * @returns {string}
 */
function lighten(hex, factor = 0.5) {
    const num = parseInt(hex.slice(1), 16);
    let r = (num >> 16) & 0xff;
    let g = (num >> 8) & 0xff;
    let b = num & 0xff;
    r = Math.round(r + (255 - r) * factor);
    g = Math.round(g + (255 - g) * factor);
    b = Math.round(b + (255 - b) * factor);
    return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Build the dynamic product checkboxes (NRPS/PKS/RiPP/...) in #typeFilterGroup.
 * @param {string[]} types
 */
function renderTypeFilters(types) {
    const group = document.getElementById('typeFilterGroup');
    if (!group) return;

    // Preserve the legend if present
    const legend = group.querySelector('legend');
    group.innerHTML = '';
    if (legend) group.appendChild(legend);

    const frag = document.createDocumentFragment();
    types.forEach(type => {
        const id = `type_${type.replace(/\W+/g, '_')}`;
        const label = document.createElement('label');
        const input = document.createElement('input');
        input.type = 'checkbox';
        input.id = id;
        input.dataset.type = type;
        input.checked = true; // default: include
        label.appendChild(input);
        label.append(` ${type}`);
        frag.appendChild(label);
    });
    group.appendChild(frag);
}

/**
 * Build the dynamic completeness checkboxes (Complete/Incomplete/Unknown completeness)
 * based only on what appears in the table.
 * @param {string[]} statuses
 */
function renderCompletenessFilters(statuses) {
    const group = document.getElementById('completenessFilterGroup');
    if (!group) return;

    // Preserve the <legend>
    const legend = group.querySelector('legend');
    group.innerHTML = '';
    if (legend) group.appendChild(legend);

    const frag = document.createDocumentFragment();
    statuses.forEach(status => {
        const id = 'status_' + status.replace(/\s+/g, '_');
        const label = document.createElement('label');
        const input = document.createElement('input');
        input.type = 'checkbox';
        input.id = id;
        input.checked = true; // default ON
        label.appendChild(input);
        label.append(` ${status}`);
        frag.appendChild(label);
    });

    group.appendChild(frag);
}

/**
 * Return which product types are currently selected in the UI.
 * Missing checkboxes default to "on".
 * @param {string[]} allTypes
 * @returns {string[]}
 */
function selectedTypesFromUI(allTypes) {
    return allTypes.filter(t => {
        const el = document.getElementById(`type_${t.replace(/\W+/g, '_')}`);
        return !el || el.checked; // include if missing or checked
    });
}

/**
 * Return completeness statuses that are currently selected in the UI.
 * @returns {('complete'|'incomplete'|'unknown completeness')[]}
 */
function selectedStatusesFromUI() {
  const on = [];
  if (document.getElementById('status_complete')?.checked)   on.push('complete');
  if (document.getElementById('status_incomplete')?.checked) on.push('incomplete');
  if (document.getElementById('status_unknown_completeness')?.checked) on.push('unknown completeness');
  // If none selected, return empty array (chart will show nothing for completeness parts)
  return on;
}

/**
 * Enable/disable controls depending on bar chart mode and "show by" choices.
 */
function syncShowByDisabled() {
    const isTotal   = document.getElementById('barModeTotal').checked;
    const showBy    = document.getElementById('barModeShowBy').checked;

    const byType = document.getElementById('byType');
    const byCompleteness = document.getElementById('byCompleteness');

    if (byType) byType.disabled = isTotal;
    if (byCompleteness) byCompleteness.disabled = isTotal;

    const byTypeOn  = byType ? byType.checked : false;
    const byCompOn  = byCompleteness ? byCompleteness.checked : false;

    const showByGroup = document.getElementById('showByGroup');
    if (showByGroup) showByGroup.disabled = isTotal;

    const typeFilterGroup = document.getElementById('typeFilterGroup');
    if (typeFilterGroup) typeFilterGroup.disabled = isTotal || !showBy || !byTypeOn;

    const completenessFilterGroup = document.getElementById('completenessFilterGroup');
    if (completenessFilterGroup) completenessFilterGroup.disabled = isTotal || !showBy || !byCompOn;
}

/** Keep an empty variable so we can update/destroy the chart cleanly. */
let bgcChart = null;

/**
 * Find a row in the data by its label (first column).
 * @param {Array[]} data
 * @param {string} label
 * @returns {Array|null}
 */
function getRowByLabel(data, label) {
    return data.find(r => r[0] === label) || null;
}

/**
 * Return all product types that exist in the table (excluding "total"/pure
 * completeness rows). Types are parsed from labels like "# BGCs (NRPS, Complete)".
 * @param {Array[]} data
 * @returns {string[]}
 */
function detectTypes(data) {
  const types = new Set();

  for (const r of data) {
    const s = String(r[0] ?? '');
    // Match labels that start with "# BGCs ( ... )"
    const m = /^#\s*BGCs?\s*\(\s*([^)]+)\s*\)\s*$/i.exec(s);
    if (!m) continue;

    // inside is: "NRPS" or "NRPS, Complete" or "Complete" etc.
    const inside = m[1].trim();

    // take only the part before the first comma (the product)
    const [typeRaw] = inside.split(',', 1);
    const type = (typeRaw || '').trim();
    if (!type) continue;

    // skip non-product tokens that appear as totals
    const lc = type.toLowerCase();
    if (lc === 'complete' || lc === 'incomplete' || lc === 'unknown completeness') continue;

    types.add(type);
  }

  // Preferred order + fallback to 'Other' if nothing detected
  const preferred = ['NRPS','PKS','RiPP','terpene','saccharide','alkaloid','hybrid'];
  const detected = Array.from(types);
  const ordered = [
    ...preferred.filter(t => detected.includes(t)),
    ...detected.filter(t => !preferred.includes(t))
  ];
  return ordered.length ? ordered : ['other'];
}

/**
 * Detect which completeness statuses exist in the table.
 * Looks for rows like "# BGCs (Complete)" or "# BGCs (Type, Unknown completeness)".
 * @param {Array[]} data
 * @returns {string[]}
 */
function detectCompleteness(data) {
    const allowed = ['complete', 'incomplete', 'unknown completeness'];
    const set = new Set();

    data.forEach(row => {
        const label = String(row[0] || '').trim();

        // total rows, e.g. "# BGCs (Unknown completeness)"
        let m = /^#\s*BGCs\s*\(\s*([^)]+)\s*\)$/i.exec(label);
        if (m) {
            const inside = m[1]
                .split(',')
                .map(s => s.trim().toLowerCase());

            inside.forEach(v => {
                if (allowed.includes(v)) {
                    set.add(v);
                }
            });
        }

        // typed rows, e.g. "# BGCs (NRPS, Unknown completeness)"
        m = /^#\s*BGCs\s*\(\s*([^,]+),\s*([^)]+)\)$/i.exec(label);
        if (m) {
            const status = m[2].trim().toLowerCase();
            if (allowed.includes(status)) {
                set.add(status);
            }
        }
    });

    return Array.from(set);
}

// Which metric labels exist for each running mode
const METRIC_TABS_BY_MODE = {
    compare_to_reference: ['bgcs', 'fully', 'partial', 'missed'],
    compare_tools:        ['bgcs', 'unique', 'pyplots'],
    compare_samples:      ['bgcs']  // Overview only
};

// Which metric are we plotting (# BGCs, Fully recovered, etc.)
const METRIC_BASE = {
    bgcs:   '# BGCs',
    fully:  '# full ref. BGCs',
    partial:'# partial ref. BGCs',
    missed: '# missed ref. BGCs',
    unique: '# unique BGCs',
    // pyplots: no metricBase → handled specially
};

let currentMetricKey = 'bgcs';

function metricLabel(base, inside) {
    // e.g. base="# BGCs", inside="Total"
    return inside ? `${base} (${inside})` : base;
}


function renderExternalLegend(chart) {
    const legendContainer = document.getElementById('bgcLegend');
    if (!legendContainer) return;

    legendContainer.innerHTML = '';

    chart.data.datasets.forEach(ds => {
        const values = Array.isArray(ds.data) ? ds.data : [];
        const hasNonZero = values.some(v => v !== 0);

        // Skip legend entries for datasets that are entirely zero
        if (!hasNonZero) return;

        const item = document.createElement('div');
        item.style.display = 'flex';
        item.style.alignItems = 'center';
        item.style.gap = '6px';

        const swatch = document.createElement('span');
        swatch.className = 'plot-legend-swatch';
        swatch.style.backgroundColor = ds.backgroundColor;

        const label = document.createElement('span');
        label.textContent = ds.label;

        item.appendChild(swatch);
        item.appendChild(label);
        legendContainer.appendChild(item);
    });
}

/**
 * Build or rebuild the bar plot based on the current UI state.
 * @param {Array[]} data
 */
function buildBarPlotDynamic(data) {
    // If the current tab does not use a bar chart (e.g. "pyplots"), do nothing
    if (!METRIC_BASE[currentMetricKey]) {
        return;
    }

    const metricBase = METRIC_BASE[currentMetricKey];

    // In COMPARE_TO_REFERENCE, drop the reference column from the plot for all tabs except Overview (totals)
    const hideReferenceForThisMetric = isReferenceMode && currentMetricKey !== 'bgcs';
    const colStart = hideReferenceForThisMetric ? 2 : 1;

    // labels are Assembly+Tool
    const tools = data[1].slice(colStart);
    const assemblies = data[0].slice(colStart);
    const labels = tools.map((tool, i) => `${tool}\n${assemblies[i]}`);

    // read UI
    const mode = document.getElementById('barModeShowBy').checked ? 'showby' : 'total';
    const byCompleteness = document.getElementById('byCompleteness').checked;
    const byType = document.getElementById('byType').checked;

    let datasets = [];

    if (mode === 'total') {
        const row = getRowByLabel(data, metricLabel(metricBase, ''));
        if (row) {
            const counts = row.slice(colStart).map(v => parseInt(v, 10));
            datasets.push({
                label: 'Total',
                data: counts,
                backgroundColor: '#322b7a'
            });
        }
    } else {
        // SHOW BY...
        if (byCompleteness && byType) {
            // stacks are (type x completeness)
            const allTypes = detectTypes(data);
            const types = selectedTypesFromUI(allTypes);
            const statuses = selectedStatusesFromUI();

            for (const type of types) {
                const baseColor = productColors[type] || productColors['Other'];
                // const rowComplete = getRowByLabel(data, `# BGCs (${type}, Complete)`);
                const rowComplete = getRowByLabel(data, metricLabel(metricBase, `${type}, complete`));
                const rowIncomplete = getRowByLabel(data, metricLabel(metricBase, `${type}, incomplete`));
                const rowUnknown    = getRowByLabel(data, metricLabel(metricBase, `${type}, unknown completeness`));

                if (statuses.includes('complete') && rowComplete) {
                    datasets.push({
                        label: `${type}, complete`,
                        data: rowComplete.slice(colStart).map(v => parseInt(v, 10)),
                        backgroundColor: baseColor
                    });
                }
                if (statuses.includes('incomplete') && rowIncomplete) {
                    datasets.push({
                        label: `${type}, incomplete`,
                        data: rowIncomplete.slice(colStart).map(v => parseInt(v, 10)),
                        backgroundColor: lighten(baseColor, 0.45) // lighter shade
                    });
                }
                if (statuses.includes('unknown completeness') && rowUnknown) {
                    datasets.push({
                        label: `${type}, unknown completeness`,
                        data: rowUnknown.slice(colStart).map(v => parseInt(v, 10)),
                        backgroundColor: lighten(baseColor, 0.75) // even lighter shade
                    });
                }
            }

        } else if (byCompleteness) {
            // stacks are complete vs incomplete (total)
            const rowComplete   = getRowByLabel(data, metricLabel(metricBase, 'complete'));
            const rowIncomplete = getRowByLabel(data, metricLabel(metricBase, 'incomplete'));
            const rowUnknown    = getRowByLabel(data, metricLabel(metricBase, 'unknown completeness'));
            const statuses = selectedStatusesFromUI();

            if (statuses.includes('complete') && rowComplete) {
                datasets.push({
                    label: 'complete',
                    data: rowComplete.slice(colStart).map(v => parseInt(v, 10)),
                    backgroundColor: '#578c18'
                });
            }
            if (statuses.includes('incomplete') && rowIncomplete) {
                datasets.push({
                    label: 'incomplete',
                    data: rowIncomplete.slice(colStart).map(v => parseInt(v, 10)),
                    backgroundColor: '#ccca3d'
                });
            }
            if (statuses.includes('unknown completeness') && rowUnknown) {
                datasets.push({
                    label: 'unknown completeness',
                    data: rowUnknown.slice(colStart).map(v => parseInt(v, 10)),
                    backgroundColor: '#999999'
                });
            }

        } else if (byType) {
            // stacks are by type (totals per type)
            const allTypes = detectTypes(data);
            const types = selectedTypesFromUI(allTypes);
            for (const type of types) {
                const row = getRowByLabel(data, metricLabel(metricBase, type));
                if (row) {
                    datasets.push({
                        label: type,
                        data: row.slice(colStart).map(v => parseInt(v, 10)),
                        backgroundColor: productColors[type] || productColors['other']
                    });
                }
            }
        } else {
            // fallback: if user selected "Show by" but no boxes,
            // just show the total to avoid an empty chart.
            const row = getRowByLabel(data, metricLabel(metricBase, ''));
            if (row) {
                const counts = row.slice(colStart).map(v => parseInt(v, 10));
                datasets.push({
                    label: 'Total',
                    data: counts,
                    backgroundColor: '#322b7a'
                });
            }
        }
    }

    // cleanup previous chart
    if (bgcChart) {
        bgcChart.destroy();  // clean up the old chart instance
        bgcChart = null;
    }

    // stacked if more than one dataset
    const stacked = datasets.length > 1;


    bgcChart = new Chart(document.getElementById('bgcBarPlot'), { // save the new chart in the same variable
        type: 'bar',
        data: { labels, datasets },
        options: {
            responsive: false,
            maintainAspectRatio: false,
            devicePixelRatio: 2,
            layout: {
                padding: { top: 10, right: 180, bottom: 10, left: 0 }
            },
            animation: {
                duration: 800,
                easing: 'easeOutQuart'
            },

            transitions: {
                // don't animate on resize at all,
                // in case something does trigger a resize
                resize: { animation: { duration: 0 } }
            },

            plugins: {
                title: {
                    display: true,
                    text: '# BGCs',
                    font: {
                        size: 16
                    },
                    padding: {
                        top: 0,
                        bottom: 10
                    }
                },
                legend: {
                    display: false
                },
                tooltip: {
                    animation: { duration: 0 },
                    callbacks: {
                        label: ctx => `${ctx.dataset.label}: ${ctx.raw}`
                    }
                }
            },

            scales: {
                x: {
                    stacked,
                    ticks: {
                        autoSkip: true,
                        autoSkipPadding: 10,
                        maxRotation: 45,
                        minRotation: 45
                    }
                },
                y: { beginAtZero: true, stacked }
            }
        }
    });
    renderExternalLegend(bgcChart);
}



/* ---------------------------------------------------------------------------
 * VENN DIAGRAM BUILDER
 * ------------------------------------------------------------------------- */
// Normalize tool names so they display consistently in the Venn diagram
function normalizeToolLabel(name) {
    if (!name) return '';
    const s = String(name).trim();
    const low = s.toLowerCase();
    if (low.startsWith('deepbgc'))   return 'DeepBGC';
    if (low.startsWith('antismash')) return 'antiSMASH';
    if (low.startsWith('gecco'))     return 'GECCO';
    return s;
}

// Draw a 2-set Venn diagram inside the given <svg> element.
// The diagram uses:
//  - left circle  = unique + non_unique from A→B
//  - right circle = unique + non_unique from B→A
//  - intersection = “A_non_unique / B_non_unique”
function drawVenn(svg, toolA, toolB, pairwiseByTool, threshold) {
    const pair = pairwiseByTool || {};

    const fromA = (pair[toolA] && pair[toolA][toolB]) || {};
    const fromB = (pair[toolB] && pair[toolB][toolA]) || {};

    const leftUnique  = Number(fromA.unique     || 0);
    const leftNonUni  = Number(fromA.non_unique || 0);
    const rightUnique = Number(fromB.unique     || 0);
    const rightNonUni = Number(fromB.non_unique || 0);

    const labelA = normalizeToolLabel(toolA);
    const labelB = normalizeToolLabel(toolB);

    // Clear the SVG before drawing a new diagram
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    svg.setAttribute('viewBox', '0 0 290 250');

    const NS = 'http://www.w3.org/2000/svg';

    const makeCircle = (cx, cy, r, fill) => {
        const c = document.createElementNS(NS, 'circle');
        c.setAttribute('cx', cx);
        c.setAttribute('cy', cy);
        c.setAttribute('r', r);
        c.setAttribute('fill', fill);
        c.setAttribute('fill-opacity', '0.7');
        c.setAttribute('stroke', '#333333');
        c.setAttribute('stroke-width', '1');
        svg.appendChild(c);
    };

    const makeText = (x, y, text, size = 16) => {
        const t = document.createElementNS(NS, 'text');
        t.setAttribute('x', x);
        t.setAttribute('y', y);
        t.setAttribute('text-anchor', 'middle');
        t.setAttribute('dominant-baseline', 'middle');
        t.setAttribute('font-size', size);
        t.setAttribute('font-family', 'Arial, sans-serif');
        t.textContent = text;
        svg.appendChild(t);
    };

    // circles
    makeCircle(115, 110, 75, '#2e808f');  // left
    makeCircle(175, 110, 75, '#FFBC42');  // right

    // counts
    makeText(85, 110,  String(leftUnique));
    makeText(205, 110, String(rightUnique));
    makeText(145, 110, `${leftNonUni} | ${rightNonUni}`);

    // labels
    makeText(75, 200, labelA, 9);
    makeText(215, 200, labelB, 9);

    // title
    const wrapper = svg.closest('.venn-wrapper');
    if (wrapper) {
        const titleEl = wrapper.querySelector('.venn-title');
        if (titleEl) {
            const thrPct = (typeof threshold === 'number') ? Math.round(threshold * 100) : null;
            titleEl.textContent =
                `Overlap between ${labelA} and ${labelB}` +
                (thrPct != null ? ` (overlap threshold = ${thrPct}%)` : '');
        }
    }
}

// Create the full Venn panel (checkboxes + title + SVG).
function initVennPanel(panel, metadata) {
    const pairwise = metadata && metadata.pairwise_by_run;
    if (!pairwise) {
        panel.textContent = 'No overlap information available.';
        return;
    }

    const tools = Object.keys(pairwise);
    if (!tools.length) {
        panel.textContent = 'No overlap information available.';
        return;
    }

    panel.innerHTML = '';

    // Container for all Venn-related UI
    const wrapper = document.createElement('div');
    wrapper.className = 'venn-wrapper';

    // Checkbox controls
    const controls = document.createElement('fieldset');
    controls.className = 'venn-controls';

    const label = document.createElement('legend');
    label.textContent = 'Genome mining tools (max. 2)';
    controls.appendChild(label);

    const boxContainer = document.createElement('span');
    boxContainer.className = 'venn-checkboxes';
    controls.appendChild(boxContainer);

    // Title above the SVG (also updated by drawVenn)
    const title = document.createElement('div');
    title.className = 'venn-title';
    title.textContent = 'Select two tools to see overlap.';

    // The SVG that will hold the circles/text
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.classList.add('venn-svg');

    // Track which tools are selected, and the current filename for download
    const selected = [];
    let currentFilename = 'bgc-quast_venn.png';

    // Download button for the Venn plot
    const downloadBtn = document.createElement('button');
    downloadBtn.id = 'downloadVennPlot';
    downloadBtn.textContent = 'Download plot';
    downloadBtn.addEventListener('click', () => {
        // Do nothing if no diagram yet (no tools selected or not 2 selected)
        if (!svg.firstChild) return;
        exportSvg(svg, currentFilename);
    });

    // Assemble UI
    // Left column: the SVG
    const mainRow = document.createElement('div');
    mainRow.className = 'venn-main-row';

    const leftCol = document.createElement('div');
    leftCol.className = 'venn-left';
    leftCol.appendChild(svg);

    const rightCol = document.createElement('div');
    rightCol.className = 'venn-right';
    rightCol.appendChild(controls);

    mainRow.appendChild(leftCol);
    mainRow.appendChild(rightCol);

    wrapper.appendChild(title);
    wrapper.appendChild(mainRow);

    // Insert download button BELOW the legend inside the right column
    rightCol.appendChild(downloadBtn);

    // Finally attach wrapper to panel
    panel.appendChild(wrapper);

    // Create one checkbox per tool
    tools.forEach(tool => {
        const lbl = document.createElement('label');
        const cb  = document.createElement('input');
        cb.type   = 'checkbox';
        cb.value  = tool;

        const txt = document.createElement('span');
        txt.textContent = normalizeToolLabel(tool);

        lbl.appendChild(cb);
        lbl.appendChild(txt);
        boxContainer.appendChild(lbl);

        cb.addEventListener('change', () => {
            if (cb.checked) {
                if (selected.length >= 2) {
                    cb.checked = false;
                    return;
                }
                selected.push(tool);
            } else {
                const idx = selected.indexOf(tool);
                if (idx !== -1) selected.splice(idx, 1);
            }

            const all = boxContainer.querySelectorAll('input[type="checkbox"]');
            all.forEach(other => {
                if (!other.checked) {
                    other.disabled = (selected.length >= 2);
                } else {
                    other.disabled = false;
                }
            });

            if (selected.length === 2) {
                drawVenn(
                    svg,
                    selected[0],
                    selected[1],
                    pairwise,
                    metadata && metadata.compare_tools_overlap_threshold
                );

                const thr = metadata && metadata.compare_tools_overlap_threshold;
                const thrPct = (typeof thr === 'number') ? Math.round(thr * 100) : null;

                const safeA = normalizeToolLabel(selected[0]).replace(/\s+/g, '');
                const safeB = normalizeToolLabel(selected[1]).replace(/\s+/g, '');

                currentFilename = (thrPct != null)
                    ? `venn_${safeA}_vs_${safeB}_${thrPct}.png`
                    : `venn_${safeA}_vs_${safeB}.png`;

            } else {
                while (svg.firstChild) svg.removeChild(svg.firstChild);
                title.textContent = 'Select two tools to see overlap.';
                currentFilename = 'bgc-quast_venn.png';
            }
        });
    });
}
/* ---------------------------------------------------------------------------
 * GCF VENN DIAGRAM (BiG-SCAPE, compare-samples)
 *
 * Separate from drawVenn above, which is a 2-set, tool-labelled, count-driven
 * compare-tools widget. This one takes real set membership for 2 or 3 samples.
 * ------------------------------------------------------------------------- */

const GCF_VENN_COLORS = ['#2e808f', '#FFBC42', '#7B6CB0'];

// Split family IDs into the exclusive / pairwise / triple regions of a Venn.
function gcfVennRegions(sets, labels) {
    const asSet = labels.map(l => new Set(sets[l] || []));
    const inOnly = (i) => [...asSet[i]].filter(
        f => asSet.every((s, j) => j === i || !s.has(f))
    ).length;
    const inBoth = (i, j) => [...asSet[i]].filter(
        f => asSet[j].has(f) && asSet.every((s, k) => (k === i || k === j) || !s.has(f))
    ).length;
    const inAll = () => [...asSet[0]].filter(f => asSet.every(s => s.has(f))).length;

    return { asSet, inOnly, inBoth, inAll };
}

// Draw a 2- or 3-set Venn of gene cluster families into the given <svg>.
function drawVennGcf(svg, sets, labels) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);

    const NS = 'http://www.w3.org/2000/svg';

    const makeCircle = (cx, cy, r, fill) => {
        const c = document.createElementNS(NS, 'circle');
        c.setAttribute('cx', cx);
        c.setAttribute('cy', cy);
        c.setAttribute('r', r);
        c.setAttribute('fill', fill);
        c.setAttribute('fill-opacity', '0.55');
        c.setAttribute('stroke', '#333333');
        c.setAttribute('stroke-width', '1');
        svg.appendChild(c);
    };

    const makeText = (x, y, text, size = 15, weight = 'normal') => {
        const t = document.createElementNS(NS, 'text');
        t.setAttribute('x', x);
        t.setAttribute('y', y);
        t.setAttribute('text-anchor', 'middle');
        t.setAttribute('dominant-baseline', 'middle');
        t.setAttribute('font-size', size);
        t.setAttribute('font-weight', weight);
        t.setAttribute('font-family', 'Arial, sans-serif');
        t.textContent = text;
        svg.appendChild(t);
    };

    if (labels.length < 2 || labels.length > 3) {
        svg.setAttribute('viewBox', '0 0 300 120');
        makeText(150, 60, `${labels.length} samples — Venn shows 2 or 3`, 13);
        return;
    }

    const { inOnly, inBoth, inAll } = gcfVennRegions(sets, labels);

    if (labels.length === 2) {
        svg.setAttribute('viewBox', '0 0 300 240');
        makeCircle(120, 110, 75, GCF_VENN_COLORS[0]);
        makeCircle(180, 110, 75, GCF_VENN_COLORS[1]);

        makeText(85, 110, String(inOnly(0)));
        makeText(215, 110, String(inOnly(1)));
        makeText(150, 110, String(inBoth(0, 1)), 15, 'bold');

        makeText(85, 215, labels[0], 11);
        makeText(215, 215, labels[1], 11);
        return;
    }

    // Three circles, centres on an equilateral triangle.
    svg.setAttribute('viewBox', '0 0 300 290');
    makeCircle(115, 110, 72, GCF_VENN_COLORS[0]);
    makeCircle(185, 110, 72, GCF_VENN_COLORS[1]);
    makeCircle(150, 175, 72, GCF_VENN_COLORS[2]);

    // exclusive
    makeText(75, 85, String(inOnly(0)));
    makeText(225, 85, String(inOnly(1)));
    makeText(150, 225, String(inOnly(2)));
    // pairwise
    makeText(150, 78, String(inBoth(0, 1)));
    makeText(103, 165, String(inBoth(0, 2)));
    makeText(197, 165, String(inBoth(1, 2)));
    // all three
    makeText(150, 137, String(inAll()), 16, 'bold');

    // labels outside the circles
    makeText(62, 25, labels[0], 11);
    makeText(238, 25, labels[1], 11);
    makeText(150, 275, labels[2], 11);
}

// Small read-only table of the GCF metrics for the selected cutoff.
function buildGcfSummaryTable(bigscape, cutoff) {
    const table = document.createElement('table');
    table.className = 'gcf-summary-table';

    const head = document.createElement('tr');
    head.appendChild(document.createElement('th'));
    (bigscape.columns || []).forEach(label => {
        const th = document.createElement('th');
        th.textContent = label;
        head.appendChild(th);
    });
    table.appendChild(head);

    (bigscape.rows[cutoff] || []).forEach(row => {
        const tr = document.createElement('tr');
        row.forEach((cell, idx) => {
            const el = document.createElement(idx === 0 ? 'th' : 'td');
            el.textContent = cell;
            tr.appendChild(el);
        });
        table.appendChild(tr);
    });

    return table;
}

// Build the whole BiG-SCAPE panel: cutoff dropdown, Venn, summary table, link out.
function initGcfPanel(panel, bigscape) {
    if (!bigscape || !bigscape.cutoffs || !bigscape.cutoffs.length) {
        panel.textContent = 'No BiG-SCAPE results available.';
        return;
    }

    panel.innerHTML = '';

    const labels = bigscape.columns || [];
    let cutoff = bigscape.default_cutoff || bigscape.cutoffs[0];

    const wrapper = document.createElement('div');
    wrapper.className = 'venn-wrapper';

    const title = document.createElement('div');
    title.className = 'venn-title';

    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.classList.add('venn-svg');

    // Controls: one dropdown, since a family only means something at one cutoff.
    const controls = document.createElement('fieldset');
    controls.className = 'venn-controls';

    const legend = document.createElement('legend');
    legend.textContent = 'GCF distance cutoff';
    controls.appendChild(legend);

    const select = document.createElement('select');
    select.className = 'gcf-cutoff-select';
    bigscape.cutoffs.forEach(c => {
        const opt = document.createElement('option');
        opt.value = c;
        opt.textContent = c;
        if (c === cutoff) opt.selected = true;
        select.appendChild(opt);
    });
    controls.appendChild(select);

    const tableHolder = document.createElement('div');
    tableHolder.className = 'gcf-summary-holder';

    const note = document.createElement('div');
    note.className = 'gcf-note';

    const link = document.createElement('a');
    link.className = 'gcf-report-link';
    link.href = bigscape.report_url || '#';
    link.textContent = 'Open full BiG-SCAPE report';
    link.target = '_blank';
    link.rel = 'noopener';

    const downloadBtn = document.createElement('button');
    downloadBtn.id = 'downloadGcfVenn';
    downloadBtn.textContent = 'Download plot';
    downloadBtn.addEventListener('click', () => {
        if (!svg.firstChild) return;
        exportSvg(svg, `gcf_venn_c${cutoff}.png`);
    });

    const mainRow = document.createElement('div');
    mainRow.className = 'venn-main-row';

    const leftCol = document.createElement('div');
    leftCol.className = 'venn-left';
    leftCol.appendChild(svg);

    const rightCol = document.createElement('div');
    rightCol.className = 'venn-right';
    rightCol.appendChild(controls);
    rightCol.appendChild(tableHolder);
    rightCol.appendChild(note);
    rightCol.appendChild(link);
    rightCol.appendChild(downloadBtn);

    mainRow.appendChild(leftCol);
    mainRow.appendChild(rightCol);
    wrapper.appendChild(title);
    wrapper.appendChild(mainRow);
    panel.appendChild(wrapper);

    const render = () => {
        title.textContent =
            `Gene cluster families shared between samples (cutoff ${cutoff})`;
        drawVennGcf(svg, (bigscape.sets && bigscape.sets[cutoff]) || {}, labels);

        tableHolder.innerHTML = '';
        tableHolder.appendChild(buildGcfSummaryTable(bigscape, cutoff));

        note.textContent = (cutoff === bigscape.default_cutoff)
            ? 'The report table above uses this cutoff.'
            : `The report table above still uses cutoff ${bigscape.default_cutoff}. `
              + 'Re-run bgc-quast with --bigscape-cutoff to change it.';
    };

    select.addEventListener('change', () => {
        cutoff = select.value;
        render();
    });

    render();
}


/* ---------------------------------------------------------------------------
 * BARPLOT DOWNLOADING
 * ------------------------------------------------------------------------- */
function exportPlotWithLegend() {
    const chartCanvas = document.getElementById('bgcBarPlot');
    const legendContainer = document.getElementById('bgcLegend');

    if (!chartCanvas || !legendContainer) return;
    if (!bgcChart) return;

    const legendItems = Array.from(legendContainer.querySelectorAll('div'));
    const hasLegend = legendItems.length > 0;

    // Base dimensions from the existing chart canvas
    const chartWidth = 1200;   // actual pixel width
    const chartHeight = chartCanvas.height; // actual pixel height

    // Simple layout constants for the legend area
    const legendPaddingX = 0;
    const legendPaddingY = 150;
    const legendLineHeight = 50;
    const legendSwatchSize = 14;
    const legendWidth = 650

    // Create an off-screen canvas
    const exportCanvas = document.createElement('canvas');
    exportCanvas.width = chartWidth + (hasLegend ? legendWidth + legendPaddingX : 0);
    exportCanvas.height = chartHeight;

    const ctx = exportCanvas.getContext('2d');

    // White background
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, exportCanvas.width, exportCanvas.height);

    // Draw the barplot
    ctx.drawImage(chartCanvas, 0, 0);

    // Draw the legend on the right
    if (hasLegend) {
        let x0 = chartWidth + legendPaddingX;
        let y0 = legendPaddingY;

        ctx.font = '32px Arial';
        ctx.textBaseline = 'middle';

        legendItems.forEach((item, idx) => {
            const swatch = item.querySelector('.plot-legend-swatch');
            if (!swatch) return;

            const color = window.getComputedStyle(swatch).backgroundColor;
            const labelText = item.textContent.trim();
            const y = y0 + idx * legendLineHeight;

            // Swatch
            ctx.fillStyle = color;
            ctx.fillRect(x0, y - legendSwatchSize / 2, legendSwatchSize, legendSwatchSize);
            ctx.strokeStyle = '#696969';
            ctx.strokeRect(x0, y - legendSwatchSize / 2, legendSwatchSize, legendSwatchSize);

            // Text
            ctx.fillStyle = '#000000';
            ctx.fillText(labelText, x0 + legendSwatchSize + 6, y);
        });
    }

    // Trigger browser download
    const svgData = `<svg xmlns="http://www.w3.org/2000/svg" width="${exportCanvas.width}" height="${exportCanvas.height}"> <image href="${exportCanvas.toDataURL('image/png')}" width="100%" height="100%"/></svg>`;

    const blob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);

    const link = document.createElement('a');
    link.href = url;
    link.download = 'bgc-quast_barplot.svg';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

/* ---------------------------------------------------------------------------
 * VENN DIAGRAM DOWNLOADING
 * ------------------------------------------------------------------------- */
// Export a single SVG element (e.g. the Venn diagram) as an SVG file
function exportSvg(svgElement, filename) {
    if (!svgElement) return;

    const serializer = new XMLSerializer();
    let source = serializer.serializeToString(svgElement);

    // Add XML header
    source = '<?xml version="1.0" standalone="no"?>\r\n' + source;

    const blob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);

    const link = document.createElement('a');
    link.href = url;
    link.download = filename.replace('.png', '.svg');
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    URL.revokeObjectURL(url);
}




/* ---------------------------------------------------------------------------
 * PAGE INITIALIZATION & EVENT WIRING
 * ------------------------------------------------------------------------- */
document.addEventListener('DOMContentLoaded', () => {
    if (typeof Chart !== 'undefined') {
        Chart.defaults.font.size = 13;
        Chart.defaults.font.family = "'Arial', sans-serif";
    }
    renderRunSummary(reportMetadata);
    buildTable(reportData);
    renderTypeFilters(detectTypes(reportData));
    renderCompletenessFilters(detectCompleteness(reportData));

    buildBarPlotDynamic(reportData);

    // Extended report toggle visibility (show/hide extra rows)
    const productBtn = document.getElementById('toggleProductBtn');
    if (productBtn) {
        productBtn.addEventListener('click', () => {
            showProductRows = !showProductRows;
            applyExtendedState();
        });
    }

    const completenessBtn = document.getElementById('toggleCompletenessBtn');
    if (completenessBtn) {
        completenessBtn.addEventListener('click', () => {
            showCompletenessRows = !showCompletenessRows;
            applyExtendedState();
        });
    }

    applyExtendedState()

    // Heatmap toggle
    const heatmapToggle = document.getElementById('heatmapToggle');
    heatmapToggle.addEventListener('change', () => {
        const show = heatmapToggle.checked;
        allNumericCells.forEach(({ cells, values }) => {
            if (show) {
                heatMapOneRow(cells, values);
            } else {
                cells.forEach(cell => {
                    cell.style.backgroundColor = '';
                    cell.style.color = '';
                });
            }
        });
    });

    // Include-reference toggle
    const includeRefToggle = document.getElementById('heatmapIncludeRef');
    if (includeRefToggle) {
        const includeRefLabel = includeRefToggle.closest('label');
        if (!isReferenceMode) {
            // Hide the whole label+checkbox outside reference mode
            if (includeRefLabel) includeRefLabel.style.display = 'none';
        } else {
            includeRefToggle.addEventListener('change', () => {
                includeReferenceInHeatmap = includeRefToggle.checked;
                const container = document.getElementById('reportTableContainer');
                container.innerHTML = '';
                buildTable(reportData);
                applyExtendedState(); // <-- keep extended state
            });
        }
    }


    // Metric tabs: which ones are visible depends on running mode
    const allMetricTabs = Array.from(document.querySelectorAll('.metric-tab'));
    const mode = (typeof reportMode === 'string') ? reportMode : 'compare_samples';
    const bigscapeData = reportMetadata && reportMetadata.bigscape;
    const allowedKeys = (METRIC_TABS_BY_MODE[mode] || ['bgcs']).slice();

    // compare-samples only gets the plots tab when there is BiG-SCAPE data to show,
    // so a normal run looks exactly as it did before.
    if (mode === 'compare_samples' && bigscapeData && !allowedKeys.includes('pyplots')) {
        allowedKeys.push('pyplots');
        const pyBtn = allMetricTabs.find(b => b.dataset.metric === 'pyplots');
        if (pyBtn) pyBtn.textContent = 'GCF overlap';
    }

    // Populate pythonPlotsPanel
    const pyPanel = document.getElementById('pythonPlotsPanel');
    if (mode === 'compare_tools') {
        if (pyPanel && reportMetadata) {
            initVennPanel(pyPanel, reportMetadata);
        }
    } else if (mode === 'compare_samples' && bigscapeData) {
        if (pyPanel) {
            initGcfPanel(pyPanel, bigscapeData);
        }
    }

    // Show/hide buttons according to mode
    allMetricTabs.forEach(btn => {
        const key = btn.dataset.metric;
        const visible = allowedKeys.includes(key);
        btn.style.display = visible ? '' : 'none';
    });

    // The default active tab is always "Overview"
    currentMetricKey = 'bgcs';
    const overviewBtn = allMetricTabs.find(b => b.dataset.metric === 'bgcs');
    if (overviewBtn) {
        allMetricTabs.forEach(b => b.classList.toggle('active', b === overviewBtn));
    }

    // Click handling (bar chart vs python plots)
    allMetricTabs.forEach(btn => {
        btn.addEventListener('click', () => {
            const key = btn.dataset.metric;
            if (!allowedKeys.includes(key)) return;

            currentMetricKey = key;

            allMetricTabs.forEach(b => b.classList.toggle('active', b === btn));

            const controls    = document.querySelector('.barplot-controls');
            const plotRow     = document.querySelector('.plot-flex-row');
            const pythonPanel = document.getElementById('pythonPlotsPanel');
            const plotContainer = document.querySelector('.plot-and-controls');

            if (key === 'pyplots') {
                // Show python-generated PNGs, hide bar chart
                if (pythonPanel) pythonPanel.style.display = 'block';
                if (controls)    controls.style.display    = 'none';
                if (plotRow)    plotRow.style.display    = 'none';
                if (plotContainer) plotContainer.classList.add('pyplots-mode');

                // Optional: destroy existing chart
                if (bgcChart) {
                    bgcChart.destroy();
                    bgcChart = null;
                }
            } else {
                // Show bar chart, hide python PNGs
                if (pythonPanel) pythonPanel.style.display = 'none';
                if (controls)    controls.style.display    = '';
                if (plotRow)    plotRow.style.display    = 'flex';
                if (plotContainer) plotContainer.classList.remove('pyplots-mode');

                buildBarPlotDynamic(reportData);
            }
        });
    });


    // Bar plot controls listeners
    // "Show by" controls: radio buttons & "by ..." checkboxes
    ['barModeTotal','barModeShowBy','byCompleteness','byType'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('change', () => {
            syncShowByDisabled();
            buildBarPlotDynamic(reportData);
        });
    });

    // Product type checkbox changes
    const typeGroup = document.getElementById('typeFilterGroup');
    if (typeGroup) {
        typeGroup.addEventListener('change', (e) => {
            if (e.target && e.target.matches('input[type="checkbox"]')) {
                buildBarPlotDynamic(reportData);
            }
        });
    }

    // Completeness checkbox changes
    const compGroup = document.getElementById('completenessFilterGroup');
    if (compGroup) {
        compGroup.addEventListener('change', (e) => {
            if (e.target && e.target.matches('input[type="checkbox"]')) {
                buildBarPlotDynamic(reportData);
            }
        });
    }

    // Download button: export plot + legend as a single PNG
    const dlBtn = document.getElementById('downloadPlotWithLegend');
    if (dlBtn) {
        dlBtn.addEventListener('click', exportPlotWithLegend);
    }

    // Initialize disabled/enabled state on first load
    syncShowByDisabled();

    // Draw heatmap gradient
    drawHeatmapLegend();  // yellow to purple


});
