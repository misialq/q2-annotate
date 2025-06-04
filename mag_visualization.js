let selectedMagId = null; // Track which MAG is currently selected
let vizData = null; // Will hold the data loaded from JSON
window.gcKDEIsWeighted = true; // Default to length-weighted for KDE plot

// Utility function to parse taxonomy strings and return user-friendly names
function parseUserFriendlyTaxonomy(taxonomyString) {
  if (!taxonomyString || taxonomyString === 'N/A') {
    return 'Unknown';
  }

  // Split by semicolon and parse each level
  const levels = taxonomyString.split(';');
  const parsed = [];

  levels.forEach(level => {
    const trimmed = level.trim();
    const underscoreIndex = trimmed.indexOf('__');
    if (underscoreIndex !== -1) {
      const prefix = trimmed.substring(0, underscoreIndex);
      const name = trimmed.substring(underscoreIndex + 2).trim();
      if (name) {
        parsed.push({ prefix, name });
      }
    }
  });

  if (parsed.length === 0) {
    return 'Unknown taxonomy';
  }

  // Get the deepest (last) level
  const deepestLevel = parsed[parsed.length - 1];

  // Format based on the prefix
  if (deepestLevel.prefix === 's') {
    // For species, try to format as "Genus species" if we have genus info
    const genusLevel = parsed.find(p => p.prefix === 'g');
    if (genusLevel && deepestLevel.name.includes(' ')) {
      return deepestLevel.name; // Already formatted as "Genus species"
    } else if (genusLevel && !deepestLevel.name.includes(' ')) {
      return `${genusLevel.name} ${deepestLevel.name}`;
    } else {
      return deepestLevel.name;
    }
  } else if (deepestLevel.prefix === 'g') {
    return `${deepestLevel.name} sp.`; // "sp." indicates species not specified
  } else if (deepestLevel.prefix === 'f') {
    return `Family ${deepestLevel.name}`;
  } else if (deepestLevel.prefix === 'o') {
    return `Order ${deepestLevel.name}`;
  } else if (deepestLevel.prefix === 'c') {
    return `Class ${deepestLevel.name}`;
  } else if (deepestLevel.prefix === 'p') {
    return `Phylum ${deepestLevel.name}`;
  } else if (deepestLevel.prefix === 'k') {
    return `Kingdom ${deepestLevel.name}`;
  } else if (deepestLevel.prefix === 'd') {
    return `Domain ${deepestLevel.name}`;
  } else {
    // For any other level (subspecies, strain, etc.), return the name as-is
    return deepestLevel.name;
  }
}

// Function to format detailed taxonomy breakdown for tooltips
function formatDetailedTaxonomy(taxonomyString) {
  if (!taxonomyString || taxonomyString === 'N/A') {
    return '<br/><strong>Taxonomy:</strong><br/>Level 1: Unknown';
  }

  // Split by semicolon and parse each level
  const levels = taxonomyString.split(';');
  const parsed = {};

  levels.forEach(level => {
    const trimmed = level.trim();
    const underscoreIndex = trimmed.indexOf('__');
    if (underscoreIndex !== -1) {
      const prefix = trimmed.substring(0, underscoreIndex);
      const name = trimmed.substring(underscoreIndex + 2).trim();
      if (name) {
        parsed[prefix] = name;
      }
    }
  });

  // Define the hierarchical order and full names
  const hierarchyOrder = [
    { key: 'd', label: 'Domain' },
    { key: 'k', label: 'Kingdom' },
    { key: 'p', label: 'Phylum' },
    { key: 'c', label: 'Class' },
    { key: 'o', label: 'Order' },
    { key: 'f', label: 'Family' },
    { key: 'g', label: 'Genus' },
    { key: 's', label: 'Species' }
  ];

  let taxonomyLines = [];
  let levelCount = 1;

  hierarchyOrder.forEach(level => {
    if (parsed[level.key]) {
      taxonomyLines.push(`Level ${levelCount}: ${parsed[level.key]}`);
      levelCount++;
    }
  });

  if (taxonomyLines.length === 0) {
    return '<br/><strong>Taxonomy:</strong><br/>Level 1: Unknown';
  }

  return '<br/><strong>Taxonomy:</strong><br/>' + taxonomyLines.join('<br/>');
}

// Function to get the deepest taxonomic level for a MAG
function getDeepestTaxonomicLevel(taxonomyString) {
  if (!taxonomyString || taxonomyString === 'N/A') {
    return null;
  }

  // Split by semicolon and parse each level
  const levels = taxonomyString.split(';');
  const parsed = {};

  levels.forEach(level => {
    const trimmed = level.trim();
    const underscoreIndex = trimmed.indexOf('__');
    if (underscoreIndex !== -1) {
      const prefix = trimmed.substring(0, underscoreIndex);
      const name = trimmed.substring(underscoreIndex + 2).trim();
      if (name) {
        parsed[prefix] = name;
      }
    }
  });

  // Check from most specific to least specific
  const hierarchyOrder = ['s', 'g', 'f', 'o', 'c', 'p', 'k', 'd'];

  for (const level of hierarchyOrder) {
    if (parsed[level]) {
      return level;
    }
  }

  return null;
}

// Function to check if a MAG meets the minimum classification level requirement
function meetsClassificationLevel(magTaxonomy, requiredLevel) {
  if (requiredLevel === 'none') {
    return true;
  }

  const deepestLevel = getDeepestTaxonomicLevel(magTaxonomy);
  if (!deepestLevel) {
    return false;
  }

  // Define hierarchy with numeric values (higher number = deeper level)
  const levelHierarchy = {
    'd': 1,  // Domain
    'k': 2,  // Kingdom
    'p': 3,  // Phylum
    'c': 4,  // Class
    'o': 5,  // Order
    'f': 6,  // Family
    'g': 7,  // Genus
    's': 8   // Species
  };

  const magLevelNumber = levelHierarchy[deepestLevel];
  const requiredLevelNumber = parseInt(requiredLevel);

  return magLevelNumber <= requiredLevelNumber;
}

const svg = d3.select('svg.chart');
const tooltip = d3.select('.tooltip');
const detailsTitle = d3.select('#details-title');
const detailsMetricsText = d3.select('#details-metrics-text');
const detailsSVG = d3.select('#details-chart');
const detailsPlaceholder = d3.select('#details-placeholder');
const warning = d3.select('#warning');
const margin = { top: 25, right: 12, bottom: 10, left: 12 }; // Increased margins further
const sampleDropdown = d3.select('#sample-select');
let activeQualityFilters = new Set(); // For multi-select quality filtering

// Contamination Slider
const contThresholdSlider = d3.select('#cont-threshold');
const contThresholdValueDisplay = d3.select('#cont-threshold-value-display');

// Completeness Slider
const compThresholdSlider = d3.select('#comp-threshold');
const compThresholdValueDisplay = d3.select('#comp-threshold-value-display');

function updateContaminationValueDisplay() {
  const value = contThresholdSlider.property('value');
  contThresholdValueDisplay.text(`${value}%`);
}

function updateCompletenessValueDisplay() {
  const value = compThresholdSlider.property('value');
  compThresholdValueDisplay.text(`${value}%`);
}

// Utility function to get ancestor taxonomy at a specific level
function getAncestorTaxonomy(fullTaxonomyString, level) {
  if (!fullTaxonomyString || level === "none" || !level) {
    return fullTaxonomyString;
  }
  const levelInt = parseInt(level, 10);
  if (isNaN(levelInt) || levelInt <= 0) {
    return fullTaxonomyString;
  }
  const parts = fullTaxonomyString.split(';');
  if (levelInt >= parts.length) {
    return fullTaxonomyString;
  }
  return parts.slice(0, levelInt).join(';');
}

// Per-MAG monochromatic color scale using accent color variations
function makeColorScale(taxa) {
  if (!taxa || taxa.length === 0) {
    return d3.scaleOrdinal().range(["#ccc"]);
  }

  // Get the base accent color from CSS
  const baseColor = getComputedStyle(document.documentElement).getPropertyValue('--accent-primary').trim();

  // Convert hex to HSL for easier manipulation
  function hexToHsl(hex) {
    let r, g, b;
    if (hex.length === 4) {
      r = parseInt(hex[1] + hex[1], 16) / 255;
      g = parseInt(hex[2] + hex[2], 16) / 255;
      b = parseInt(hex[3] + hex[3], 16) / 255;
    } else {
      r = parseInt(hex.slice(1, 3), 16) / 255;
      g = parseInt(hex.slice(3, 5), 16) / 255;
      b = parseInt(hex.slice(5, 7), 16) / 255;
    }

    const max = Math.max(r, g, b);
    const min = Math.min(r, g, b);
    let h, s, l = (max + min) / 2;

    if (max === min) {
      h = s = 0;
    } else {
      const d = max - min;
      s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
      switch (max) {
        case r: h = (g - b) / d + (g < b ? 6 : 0); break;
        case g: h = (b - r) / d + 2; break;
        case b: h = (r - g) / d + 4; break;
      }
      h /= 6;
    }

    return [h * 360, s * 100, l * 100];
  }

  const [baseHue, baseSaturation, baseLightness] = hexToHsl(baseColor);
  const numTaxa = taxa.length;

  const colors = taxa.map((taxon, index) => {
    if (numTaxa === 1) {
      return `hsl(${baseHue}, ${baseSaturation}%, ${baseLightness}%)`;
    }

    // Create variations by adjusting saturation and lightness
    // Most important (index 0) gets darkest, highest saturation
    // Less important gets lighter, less saturated

    const progress = index / (numTaxa - 1);

    // Saturation: 90% to 40% (high to low importance)
    const saturation = 90 - (progress * 50);

    // Lightness: 25% to 70% (dark to light, dark = more important)
    const lightness = 25 + (progress * 45);

    return `hsl(${baseHue}, ${saturation}%, ${lightness}%)`;
  });

  return d3.scaleOrdinal().domain(taxa).range(colors);
}

function renderChart() {
  if (!vizData) {
    console.error("Visualization data not loaded yet!");
    return;
  }
  const selectedSample = sampleDropdown.property('value');
  const rawMags = vizData[selectedSample] || [];

  // Apply quality filters first
  let filteredByQualityMags = rawMags;
  if (activeQualityFilters.size > 0) {
    filteredByQualityMags = rawMags.filter(mag => activeQualityFilters.has(mag.quality));
  }

  // Apply taxonomic classification level filter
  const selectedClassificationLevel = d3.select('#classification-level-select').property('value');
  let filteredByClassificationMags = filteredByQualityMags.filter(mag =>
    meetsClassificationLevel(mag.assigned_taxonomy, selectedClassificationLevel)
  );

  if (!filteredByClassificationMags.length && filteredByQualityMags.length > 0) {
    svg.selectAll('*').remove();
    warning.text('No MAGs meet the selected taxonomic classification level requirement.').style('display', 'block');
    svg.attr('viewBox', `0 0 ${svg.node().clientWidth} ${margin.top + margin.bottom + 40}`)
       .attr('height', margin.top + margin.bottom + 40);
    updateDetailsView(null, 1); // Clear details view
    return;
  } else if (!filteredByClassificationMags.length && rawMags.length > 0 && activeQualityFilters.size > 0) {
    svg.selectAll('*').remove();
    warning.text('No MAGs match the selected quality filter(s) for this sample.').style('display', 'block');
    svg.attr('viewBox', `0 0 ${svg.node().clientWidth} ${margin.top + margin.bottom + 40}`)
       .attr('height', margin.top + margin.bottom + 40);
    updateDetailsView(null, 1); // Clear details view
    return;
  } else if (!filteredByClassificationMags.length && rawMags.length === 0) {
    svg.selectAll('*').remove();
    warning.text('No MAGs were recovered in this sample.').style('display', 'block');
    svg.attr('viewBox', `0 0 ${svg.node().clientWidth} ${margin.top + margin.bottom + 40}`)
    .attr('height', margin.top + margin.bottom + 40);
    updateDetailsView(null, 1); // Clear details view
    return;
  }
  warning.style('display', 'none');

  const sortKey = d3.select('#sort-select').property('value');
  const sortOrder = d3.select('#order-select').property('value');
  const contaminationFilter = +contThresholdSlider.property('value');
  const completenessFilter = +compThresholdSlider.property('value');
  const selectedCollapseLevel = d3.select('#collapse-level-select').property('value'); // Get collapse level

  // Now apply contamination and completeness filters to the already filtered MAGs
  let mags = filteredByClassificationMags.filter(d => 
    d.contamination <= contaminationFilter &&
    d.completeness >= completenessFilter // Changed from <= to >=
  );
  
  if (mags.length === 0 && filteredByClassificationMags.length > 0) {
    svg.selectAll('*').remove();
    warning.text('No MAGs pass the current contamination/min completeness thresholds with the selected filters.').style('display', 'block'); // Updated warning message
    svg.attr('viewBox', `0 0 ${svg.node().clientWidth} ${margin.top + margin.bottom + 40}`)
    .attr('height', margin.top + margin.bottom + 40);
    updateDetailsView(null, 1); // Clear details view
    return;
  } else if (mags.length === 0) {
    // Handle other empty cases
    updateDetailsView(null, 1); // Clear details view
    return;
  }

  warning.style('display', 'none');

  mags.sort((a, b) => {
    const valA = a[sortKey];
    const valB = b[sortKey];
    const cmp = (valA < valB) ? -1 : (valA > valB) ? 1 : 0;
    return sortOrder === 'ascending' ? cmp : -cmp;
  });

  const svgWidth = svg.node().clientWidth;
  const dataAreaWidth = svgWidth - margin.left - margin.right;

  const rowVisualHeight = 24;
  const rowGap = 8;
  const rowStep = rowVisualHeight + rowGap;

  svg.attr('viewBox', `0 0 ${svgWidth} ${mags.length * rowStep + margin.top + margin.bottom + 60}`)
     .attr('height', mags.length * rowStep + margin.top + margin.bottom + 60);
  svg.selectAll('*').remove();

  // Define layout for two plots
  const magIdWidth = 80;
  const plotGap = 20; // Gap between the two plots
  const iconColumnWidth = 30;
  const leftPlotWidth = (dataAreaWidth - magIdWidth - plotGap - iconColumnWidth) * 0.3; // 30% for quality metrics
  const rightPlotWidth = (dataAreaWidth - magIdWidth - plotGap - iconColumnWidth) * 0.7; // 70% for assembly composition

  // Calculate positions
  const leftPlotStartX = magIdWidth;
  const rightPlotStartX = leftPlotStartX + leftPlotWidth + plotGap;
  const iconXPosition = rightPlotStartX + rightPlotWidth + (iconColumnWidth / 2);

  // Left plot: contamination and completeness
  const contPortion = 0.5; // 50/50 split in left plot
  const contaminationZoneWidth = leftPlotWidth * contPortion;
  const completenessZoneWidth = leftPlotWidth * (1 - contPortion);
  const leftPlotDividerX = leftPlotStartX + contaminationZoneWidth;

  const maxContVal = d3.max(mags, d => d.contamination);
  const actualMaxContDisplay = (maxContVal === undefined || maxContVal === 0) ? 1 : maxContVal;
  const contScale = d3.scaleLinear().domain([0, actualMaxContDisplay]).range([0, contaminationZoneWidth - 10]);
  const getContBarWidth = (cVal) => {
    if (maxContVal === undefined || maxContVal === 0) return 0;
    return contScale(cVal < 0 ? 0 : cVal);
  };

  const compScale = d3.scaleLinear().domain([0, 100]).range([0, completenessZoneWidth - 10]);

  // Right plot: scale based on selected scaling method
  const scalingMethod = d3.select('#scaling-select').property('value');
  let maxRightPlotValue, rightPlotScale, rightAxisLabel;

  if (scalingMethod === 'length') {
    maxRightPlotValue = d3.max(mags, d => d.total_length);
    rightPlotScale = d3.scaleLinear().domain([0, maxRightPlotValue]).range([0, rightPlotWidth - 10]);
    rightAxisLabel = 'Assembly Length (bp)';
  } else {
    // Calculate max contig count for each MAG
    maxRightPlotValue = d3.max(mags, d => {
      const contigsMap = d.contigs || {};
      let totalContigCount = 0;
      Object.values(contigsMap).forEach(contigList => { totalContigCount += contigList.length; });
      return totalContigCount;
    });
    rightPlotScale = d3.scaleLinear().domain([0, maxRightPlotValue]).range([0, rightPlotWidth - 10]);
    rightAxisLabel = 'Contig Count';
  }

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  // Left plot axes
  g.append('g').attr('transform', `translate(0,${mags.length * rowStep})`)
    .attr('class', 'axis')
    .call(d3.axisBottom(d3.scaleLinear().domain([0,100]).range([leftPlotDividerX, leftPlotDividerX + completenessZoneWidth - 10])).ticks(Math.min(4, Math.floor(completenessZoneWidth/40))))
    .selectAll('line').attr('stroke', getComputedStyle(document.documentElement).getPropertyValue('--grid-color'));

  const contAxisScale = d3.scaleLinear().domain([0, actualMaxContDisplay]).range([leftPlotDividerX, leftPlotDividerX - getContBarWidth(actualMaxContDisplay)]);
  const contAxis = d3.axisBottom(contAxisScale);
  if (actualMaxContDisplay > 0 && actualMaxContDisplay <= 5 && Number.isInteger(actualMaxContDisplay)) {
    contAxis.ticks(actualMaxContDisplay);
  } else if (actualMaxContDisplay > 0) {
    contAxis.ticks(Math.min(4, Math.ceil(actualMaxContDisplay)));
  } else {
    contAxis.ticks(1).tickFormat(d => d === 0 ? "0" : "");
  }
  g.append('g').attr('transform', `translate(0,${mags.length * rowStep})`)
    .attr('class', 'axis')
    .call(contAxis)
    .selectAll('line').attr('stroke', getComputedStyle(document.documentElement).getPropertyValue('--grid-color'));

  // Right plot axis
  g.append('g').attr('transform', `translate(0,${mags.length * rowStep})`)
    .attr('class', 'axis')
    .call(d3.axisBottom(d3.scaleLinear().domain([0, maxRightPlotValue]).range([rightPlotStartX, rightPlotStartX + rightPlotWidth - 10]))
          .ticks(Math.min(5, Math.floor(rightPlotWidth/80)))
          .tickFormat(d3.format(".2s"))) // Use SI prefix format for large numbers
    .selectAll('line').attr('stroke', getComputedStyle(document.documentElement).getPropertyValue('--grid-color'));

  // Axis Labels
  const axisLabelYOffset = mags.length * rowStep + 35;
  g.append('text').attr('class', 'chart-axis-label').attr('x', rightPlotStartX + (rightPlotWidth / 2)).attr('y', axisLabelYOffset).attr('text-anchor', 'middle').text(rightAxisLabel);

  // Plot dividers
  g.append('line').attr('class', 'divider').attr('x1', leftPlotDividerX).attr('x2', leftPlotDividerX)
    .attr('y1', 0).attr('y2', mags.length * rowStep);

  // Rows - General Update Pattern
  const t = svg.transition().duration(750);
  const rows = g.selectAll('.row').data(mags, d => d.id);

  rows.exit().transition(t)
      .attr('transform', (d,i) => `translate(0,${i * rowStep + rowVisualHeight / 2})`)
      .style('opacity', 0).remove();

  const rowsEnter = rows.enter().append('g')
    .attr('class', 'row')
     .attr('transform', (d,i) => `translate(0,${i * rowStep - rowVisualHeight / 2})`)
     .style('opacity', 0)
     .on('click', (event, d) => {
       // Update selected MAG
       selectedMagId = d.id;

       // Update selection for all rows explicitly
       g.selectAll('.row').each(function(rowData) {
         const isSelected = rowData.id === selectedMagId;
         d3.select(this).classed('selected', isSelected);
       });
       updateDetailsView(d, actualMaxContDisplay);
     });

   // Add background rectangle for highlighting
   rowsEnter.append('rect')
     .attr('class', 'row-background')
     .attr('x', 0)
      .attr('y', 0)
     .attr('width', dataAreaWidth)
     .attr('height', rowVisualHeight)
     .attr('rx', 3);

   // Add subtle left accent bar for selection
   rowsEnter.append('rect')
     .attr('class', 'row-accent')
     .attr('x', 0)
     .attr('y', 2)
     .attr('width', 3)
     .attr('height', rowVisualHeight - 4)
     .attr('rx', 1.5);

   rowsEnter.append('text').attr('x', 5).attr('y', rowVisualHeight / 2).attr('dy', '0.35em').text(d => d.id.substring(0, 8));

   // Quality icons
   rowsEnter.each(function(d) {
    const iconFo = d3.select(this).append('foreignObject')
      .attr('x', iconXPosition - (rowVisualHeight / 2) )
      .attr('y', 0)
      .attr('width', rowVisualHeight)
      .attr('height', rowVisualHeight)
      .on('mousemove', (e) => {
        positionTooltip(e, `Quality: ${d.quality}`);
      })
      .on('mouseout', () => tooltip.style('opacity',0));

    let iconClass = 'bi '; let colorClass = '';
    if (d.quality === "High") { iconClass += 'bi-check-circle-fill'; colorClass = 'quality-high'; }
    else if (d.quality === "Medium") { iconClass += 'bi-circle-half'; colorClass = 'quality-medium'; }
    else if (d.quality === "Low") { iconClass += 'bi-exclamation-triangle-fill'; colorClass = 'quality-low'; }
    else { iconClass = ''; }
    if (iconClass) { iconFo.append('xhtml:div').style('width', '100%').style('height', '100%').style('display', 'flex').style('align-items', 'center').style('justify-content', 'center').html(`<i class="${iconClass} ${colorClass} quality-icon"></i>`); }
  });

  // Left plot bars: contamination and completeness
  rowsEnter.each(function(d) {
    const contigsMap = d.contigs || {};
    let totalContigCount = 0;
    Object.values(contigsMap).forEach(contigList => { totalContigCount += contigList.length; });

    // Contamination bar
      d3.select(this).append('rect')
      .attr('x', leftPlotDividerX - getContBarWidth(d.contamination))
      .attr('y', 0)
      .attr('width', getContBarWidth(d.contamination))
      .attr('height', rowVisualHeight)
      .attr('class', 'mag-contamination-bar').attr('rx', 6)
        .on('mousemove', e => {
        positionTooltip(e, `MAG: ${d.id}<br/>Contamination: ${d.contamination}%<br/>Completeness: ${d.completeness}%<br/>Total contigs: ${totalContigCount}`);
      }).on('mouseout', () => tooltip.style('opacity',0));

    // Completeness bar
    d3.select(this).append('rect')
      .attr('x', leftPlotDividerX + 2)
      .attr('y', 0)
      .attr('width', compScale(d.completeness))
      .attr('height', rowVisualHeight)
      .attr('fill', 'var(--success-color)')
      .attr('rx', 6)
      .on('mousemove', e => {
        positionTooltip(e, `MAG: ${d.id}<br/>Completeness: ${d.completeness}%<br/>Contamination: ${d.contamination}%<br/>Total contigs: ${totalContigCount}`);
      }).on('mouseout', () => tooltip.style('opacity',0));
  });

  // Right plot bars: contig composition scaled to assembly length
  rowsEnter.each(function(d) {
    const contigsMap = d.contigs || {};
    // const groupingTaxonomies = Object.keys(contigsMap);
    // if (groupingTaxonomies.length === 0) return; // Original check

    const scalingMethod = d3.select('#scaling-select').property('value');

    // 1. Determine original taxa order and color scale (always based on uncollapsed data)
    const originalTaxaForColoring = Object.keys(contigsMap).sort((a, b) => {
        const lengthA = d3.sum(contigsMap[a] || [], c => c.length);
        const lengthB = d3.sum(contigsMap[b] || [], c => c.length);
        return lengthB - lengthA;
    });
    const colorScaleChart = makeColorScale(originalTaxaForColoring);

    // 2. Prepare data for plotting (either collapsed or original)
    let plottedTaxaDetails = [];

    if (selectedCollapseLevel !== "none" && Object.keys(contigsMap).length > 0) {
        const aggregatedTaxa = {};
        for (const [originalGroupTaxon, contigList] of Object.entries(contigsMap)) {
            if (!contigList || contigList.length === 0) continue;
            const collapsedAncestor = getAncestorTaxonomy(originalGroupTaxon, selectedCollapseLevel);

            if (!aggregatedTaxa[collapsedAncestor]) {
                aggregatedTaxa[collapsedAncestor] = { totalLength: 0, count: 0, contributingOriginalTaxa: {} };
            }
            const currentOriginalGroupLength = d3.sum(contigList, c => c.length);
            aggregatedTaxa[collapsedAncestor].totalLength += currentOriginalGroupLength;
            aggregatedTaxa[collapsedAncestor].count += contigList.length;
            aggregatedTaxa[collapsedAncestor].contributingOriginalTaxa[originalGroupTaxon] =
                (aggregatedTaxa[collapsedAncestor].contributingOriginalTaxa[originalGroupTaxon] || 0) + currentOriginalGroupLength;
        }

        plottedTaxaDetails = Object.entries(aggregatedTaxa).map(([collapsedTaxon, data]) => {
            let dominantOriginalTaxon = "";
            let maxLen = -1;
            for (const [origTaxon, len] of Object.entries(data.contributingOriginalTaxa)) {
                if (len > maxLen) {
                    maxLen = len;
                    dominantOriginalTaxon = origTaxon;
                }
            }
            return {
                taxonDisplay: collapsedTaxon,
                colorKey: dominantOriginalTaxon || collapsedTaxon,
                totalLength: data.totalLength,
                count: data.count,
                metricValue: (scalingMethod === 'length' ? data.totalLength : data.count)
            };
        });
    } else if (Object.keys(contigsMap).length > 0) { // No collapse, but contigs exist
        plottedTaxaDetails = originalTaxaForColoring.map(groupTaxon => {
            const contigsInGroup = contigsMap[groupTaxon] || [];
            const totalLengthInGroup = d3.sum(contigsInGroup, c => c.length);
            const countInGroup = contigsInGroup.length;
            return {
                taxonDisplay: groupTaxon,
                colorKey: groupTaxon,
                totalLength: totalLengthInGroup,
                count: countInGroup,
                metricValue: (scalingMethod === 'length' ? totalLengthInGroup : countInGroup)
            };
        });
    }

    if (plottedTaxaDetails.length === 0) return; // If no data to plot for this MAG, skip

    plottedTaxaDetails.sort((a, b) => b.metricValue - a.metricValue);

    // 3. Calculate overall scale for the MAG's right plot bar
    let totalMagMetricForScaling;
    if (scalingMethod === 'length') {
        totalMagMetricForScaling = d.total_length;
    } else {
        let totalIndividualContigsInMAG = 0;
        Object.values(contigsMap).forEach(contigList => { totalIndividualContigsInMAG += (contigList || []).length; });
        totalMagMetricForScaling = totalIndividualContigsInMAG;
    }

    const availableRightWidth = rightPlotScale(totalMagMetricForScaling);
    const sumOfPlottedMetrics = d3.sum(plottedTaxaDetails, ti => ti.metricValue);

    // 4. Draw the bars
    let currentX = rightPlotStartX + 2;
    plottedTaxaDetails.forEach(taxonInfo => {
        if (taxonInfo.metricValue === 0 || sumOfPlottedMetrics === 0) return;

        const proportionOfTotalPlotted = taxonInfo.metricValue / sumOfPlottedMetrics;
        let taxonBarWidth = proportionOfTotalPlotted * availableRightWidth;

        // Ensure minimum visible width for very small segments if they are the only one or few
        // taxonBarWidth = Math.max(taxonBarWidth, 0.5);
        taxonBarWidth -= 2; // Subtract padding after potential Math.max
        if (taxonBarWidth <= 0 && proportionOfTotalPlotted > 0) taxonBarWidth = 0.5; // Miniscule bar if it has data but rounds down due to padding
        else if (taxonBarWidth <= 0) return;


        d3.select(this).append('rect').attr('x', currentX).attr('y', 0).attr('width', taxonBarWidth).attr('height', rowVisualHeight)
        .attr('fill', colorScaleChart(taxonInfo.colorKey))
        .attr('rx', 6)
        .on('mousemove', e => {
          const tooltipText = `Taxon Group: ${parseUserFriendlyTaxonomy(taxonInfo.taxonDisplay)}` +
                          `${selectedCollapseLevel !== "none" ? " (collapsed)" : ""}<br/>` +
                          `Contig Count: ${taxonInfo.count.toLocaleString()}<br/>` +
                          `Total Length in Group: ${taxonInfo.totalLength.toLocaleString()} bp<br/>` +
                          `MAG Assembly Length: ${d.total_length.toLocaleString()} bp`;
          positionTooltip(e, tooltipText);
        }).on('mouseout', () => tooltip.style('opacity',0));
      currentX += taxonBarWidth + 2;
    });
  });

  const rowsUpdate = rowsEnter.merge(rows);
  rowsUpdate.transition(t).attr('transform', (d,i) => `translate(0,${i * rowStep})`).style('opacity', 1);

  // Preserve selection state after update
  rowsUpdate.each(function(d) {
    d3.select(this).classed('selected', d.id === selectedMagId);
  });

  // Automatically update details view with the first MAG in the list
  if (mags.length > 0) {
    // If no MAG is currently selected, or the selected MAG is not in the current filtered list
    if (!selectedMagId || !mags.find(mag => mag.id === selectedMagId)) {
      selectedMagId = mags[0].id;
      updateDetailsView(mags[0], actualMaxContDisplay);
    }
    // Apply selection immediately to all rows
    rowsUpdate.each(function(d) {
      d3.select(this).classed('selected', d.id === selectedMagId);
    });
  } else {
    selectedMagId = null;
    updateDetailsView(null, actualMaxContDisplay);
  }
}

function updateDetailsView(d, actualMaxContDisplay) {
  if (!d) { // If no MAG data is provided (e.g., empty list)
    detailsTitle.text('MAG Details'); // Reset title
    detailsSVG.selectAll('*').remove(); // Clear SVG
    detailsPlaceholder.style('display', 'block'); // Show placeholder
    // Clear stored data
    window.currentMagData = null;
    window.currentMaxContDisplay = null;
    return;
  }

  // Store current data for toggle functionality
  window.currentMagData = d;
  window.currentMaxContDisplay = actualMaxContDisplay;

  detailsPlaceholder.style('display', 'none'); // Hide placeholder if we have data

  detailsTitle.text(d.id); // Just show the MAG ID

  // Populate the metrics text div with all three lines having consistent styling
  detailsMetricsText.html(
    `<table class="mag-metrics-table">
      <tr>
        <td>Assigned Taxon:</td>
        <td>${parseUserFriendlyTaxonomy(d.assigned_taxonomy)}</td>
      </tr>
      <tr>
        <td>Total Length:</td>
        <td>${d.total_length.toLocaleString()} bp</td>
      </tr>
      <tr>
        <td>N50:</td>
        <td>${d.n50.toLocaleString()} bp</td>
      </tr>
      <tr>
        <td>L50:</td>
        <td>${d.l50}</td>
      </tr>
    </table>`
  );

  detailsSVG.selectAll('*').remove();

  // Check toggle state from the top controls
  const selectedView = window.selectedDetailsView || 'circular';

  // Render initial view based on toggle state
  if (selectedView === 'bar') {
    renderBarChartView(d, actualMaxContDisplay);
  } else if (selectedView === 'gc') {
    renderMAGGCKDEPlot(d);
  } else {
    renderCircularView(d, actualMaxContDisplay);
  }
}

// Extract circular view logic to separate function
function renderCircularView(d, actualMaxContDisplay) {
  const svgNode = detailsSVG.node();
  const w = svgNode.clientWidth || 300; // Use actual width with fallback
  const h = svgNode.clientHeight || 320; // Use actual height with fallback
  const outerPlotRadius = Math.min(w, h) / 2 - 20;
  const modalBarThickness = 25;
  const coverageRingThickness = 15; // Thickness for the new coverage ring
  const ringGap = 3; // Gap between coverage ring and contig ring

  const contigsMap = d.contigs || {};

  // First, calculate total length per taxon to determine sorting order
  const taxonTotalLengths = {};
  Object.entries(contigsMap).forEach(([groupTaxon, specificContigsArray]) => {
    taxonTotalLengths[groupTaxon] = d3.sum(specificContigsArray, c => c.length);
  });

  // Create ordered list of taxa by total length (descending)
  const taxaOrderedByLength = Object.keys(taxonTotalLengths)
    .sort((a, b) => taxonTotalLengths[b] - taxonTotalLengths[a]);

  // Create shared color scale based on taxa ordered by total length
  const sharedColorScale = makeColorScale(taxaOrderedByLength);

  const allContigsFlat = Object.entries(contigsMap)
    .flatMap(([groupTaxon, specificContigsArray]) =>
        specificContigsArray.map(contig => ({
            specific_taxonomy: contig.specific_taxonomy,
            length: contig.length,
            coverage: contig.coverage || 0, // Add coverage with fallback
            group_taxonomy: groupTaxon
        }))
    )
    .filter(c => c.length > 0)
    .sort((a,b) => {
      // First sort by taxon total length (descending)
      const aTotalLength = taxonTotalLengths[a.group_taxonomy];
      const bTotalLength = taxonTotalLengths[b.group_taxonomy];
      if (aTotalLength !== bTotalLength) {
        return bTotalLength - aTotalLength;
      }
      // Then sort by individual contig length (descending) within the same group
      return b.length - a.length;
    });

  if (allContigsFlat.length === 0) {
    detailsSVG.append("text")
      .attr("x", w / 2).attr("y", h / 2)
      .attr("text-anchor", "middle").attr("dominant-baseline", "central")
      .text("No contigs with length > 0 to display.");
    return;
  }

  const totalMagContigLength = d3.sum(allContigsFlat, c => c.length);
  if (totalMagContigLength === 0) {
     detailsSVG.append("text").attr("x", w/2).attr("y", h/2).attr("text-anchor", "middle").text("Contigs have no total length.");
     return;
  }

  // Calculate max coverage for scaling the coverage ring
  const maxCoverage = d3.max(allContigsFlat, c => c.coverage);

  const completenessPercentage = d.completeness / 100;
  const totalAngleToDisplay = completenessPercentage * 2 * Math.PI;
  const arcPadding = 0.015;
  const minArcAngle = 0.02; // Minimum angle for visibility (about 1.1 degrees)

  // First, calculate what the angle would be for each contig if displayed individually
  const contigsWithAngles = allContigsFlat.map(contig => ({
    ...contig,
    wouldBeAngle: (contig.length / totalMagContigLength) * totalAngleToDisplay
  }));

  // Group contigs by taxon and separate large vs small
  const taxonGroups = d3.group(contigsWithAngles, c => c.group_taxonomy);
  const displayItems = [];

  // Process taxa in the same order (by total length)
  taxaOrderedByLength.forEach(taxon => {
    const contigs = taxonGroups.get(taxon);
    if (!contigs) return;

    const largeContigs = contigs.filter(c => c.wouldBeAngle >= minArcAngle);
    const smallContigs = contigs.filter(c => c.wouldBeAngle < minArcAngle);

    // Add large contigs individually
    largeContigs.forEach(contig => {
      displayItems.push({
        type: 'individual',
        taxon: contig.group_taxonomy,
        specific_taxonomy: contig.specific_taxonomy,
        length: contig.length,
        coverage: contig.coverage,
        count: 1
      });
    });

    // Aggregate small contigs by taxon
    if (smallContigs.length > 0) {
      const totalSmallLength = d3.sum(smallContigs, c => c.length);
      const avgCoverage = d3.mean(smallContigs, c => c.coverage);
      displayItems.push({
        type: 'aggregated',
        taxon: taxon,
        length: totalSmallLength,
        coverage: avgCoverage,
        count: smallContigs.length,
        individualLengths: smallContigs.map(c => c.length).sort((a,b) => b-a)
      });
    }
  });

  const numPads = displayItems.length > 1 ? displayItems.length : 0;
  const totalPaddingAngle = numPads * arcPadding;
  const availableAngleForData = totalAngleToDisplay - totalPaddingAngle;
  let currentAngle = 0;

  displayItems.forEach(item => {
    if (item.length === 0 || availableAngleForData <= 0) return;
    const proportionOfLength = item.length / totalMagContigLength;
    const angleForThisItem = proportionOfLength * availableAngleForData;
    const startAngle = currentAngle;
    const endAngle = currentAngle + angleForThisItem;

    // Coverage ring (outer)
    const coverageHeightProportion = (item.coverage || 0) / (maxCoverage || 1);
    const coverageRingHeight = coverageHeightProportion * coverageRingThickness;
    const coverageInnerRadius = outerPlotRadius + ringGap;
    const coverageOuterRadius = coverageInnerRadius + coverageRingHeight;

    const coverageArc = d3.arc()
      .innerRadius(coverageInnerRadius)
      .outerRadius(coverageOuterRadius)
      .startAngle(startAngle)
      .endAngle(endAngle)
      .cornerRadius(2);

    // Contig ring (inner)
    const originalArc = d3.arc().innerRadius(outerPlotRadius - modalBarThickness).outerRadius(outerPlotRadius).startAngle(startAngle).endAngle(endAngle).cornerRadius(4);
    const hoveredArc = d3.arc().innerRadius(outerPlotRadius - modalBarThickness).outerRadius(outerPlotRadius + 5).startAngle(startAngle).endAngle(endAngle).cornerRadius(4);
    const originalColor = sharedColorScale(item.taxon);

    // Draw coverage ring with grey color
    detailsSVG.append('path')
      .attr('d', coverageArc())
      .attr('transform', `translate(${w/2},${h/2})`)
      .attr('fill', '#78065C') // Deep magenta/purple
      .attr('stroke', 'white')
      .attr('stroke-width', 0.5)
      .on('mouseover', function() { 
        d3.select(this).transition().duration(150).style('fill', '#970774'); // Lighter magenta/purple on hover
        tooltip.style('opacity',1);
      })
      .on('mousemove', e => {
        let tooltipContent;
        if (item.type === 'individual') {
          tooltipContent = `<strong>Individual Contig</strong><br/>` +
                          `Primary Taxon: ${parseUserFriendlyTaxonomy(item.specific_taxonomy)}<br/>` +
                          `Group: ${parseUserFriendlyTaxonomy(item.taxon)}<br/>` +
                          `Length: ${item.length.toLocaleString()} bp<br/>` +
                          `Coverage: ${item.coverage.toFixed(1)}x` +
                          `${formatDetailedTaxonomy(item.specific_taxonomy)}`;
        } else {
          const avgLength = Math.round(item.length / item.count);
          tooltipContent = `<strong>${item.count} Small Contigs</strong><br/>` +
                          `Group: ${parseUserFriendlyTaxonomy(item.taxon)}<br/>` +
                          `Total Length: ${item.length.toLocaleString()} bp<br/>` +
                          `Average Length: ${avgLength.toLocaleString()} bp<br/>` +
                          `Average Coverage: ${item.coverage.toFixed(1)}x` +
                          `${formatDetailedTaxonomy(item.taxon)}`;
        }
        positionTooltip(e, tooltipContent);
      })
      .on('mouseout', function() { 
        d3.select(this).transition().duration(150).style('fill', '#78065C'); // Back to deep magenta/purple
        tooltip.style('opacity',0);
      });

    // Draw contig ring
    detailsSVG.append('path')
      .attr('d', originalArc())
      .attr('transform', `translate(${w/2},${h/2})`)
      .attr('fill', originalColor)
      .on('mouseover', function() {
        d3.select(this).transition().duration(150).attr('d', hoveredArc()).style('fill', d3.color(originalColor).brighter(0.7));
        tooltip.style('opacity',1);
      })
      .on('mousemove', e => {
        let tooltipContent;
        if (item.type === 'individual') {
          tooltipContent = `<strong>Individual Contig</strong><br/>` +
                          `Primary Taxon: ${parseUserFriendlyTaxonomy(item.specific_taxonomy)}<br/>` +
                          `Group: ${parseUserFriendlyTaxonomy(item.taxon)}<br/>` +
                          `Length: ${item.length.toLocaleString()} bp<br/>` +
                          `Coverage: ${item.coverage.toFixed(1)}x` +
                          `${formatDetailedTaxonomy(item.specific_taxonomy)}`;
        } else {
          const avgLength = Math.round(item.length / item.count);
          tooltipContent = `<strong>${item.count} Small Contigs</strong><br/>` +
                          `Group: ${parseUserFriendlyTaxonomy(item.taxon)}<br/>` +
                          `Total Length: ${item.length.toLocaleString()} bp<br/>` +
                          `Average Length: ${avgLength.toLocaleString()} bp<br/>` +
                          `Average Coverage: ${item.coverage.toFixed(1)}x` +
                          `${formatDetailedTaxonomy(item.taxon)}`;
        }
        positionTooltip(e, tooltipContent);
      })
      .on('mouseout', function() {
        d3.select(this).transition().duration(150).attr('d', originalArc()).style('fill', originalColor);
        tooltip.style('opacity',0);
      });
    currentAngle = endAngle + arcPadding;
  });

  const centerPadding = 10;
  const centerRadius = (outerPlotRadius - modalBarThickness) - centerPadding;
  const barWidth = 20; const barGap = 10;
  const maxBarHeight = centerRadius * 0.55;
  if (centerRadius > (barWidth * 2 + barGap)) {
    const centerChartGroup = detailsSVG.append('g')
      .attr('transform', `translate(${w/2}, ${h/2})`);
    const compHeight = (d.completeness / 100) * maxBarHeight;
    centerChartGroup.append('rect').attr('x', -barWidth - barGap / 2).attr('y', maxBarHeight - compHeight).attr('width', barWidth).attr('height', compHeight).attr('class', 'modal-completeness-bar').attr('rx', 3);
    centerChartGroup.append('text').attr('x', -barWidth / 2 - barGap / 2).attr('y', maxBarHeight - compHeight - 5).attr('text-anchor', 'middle').attr('font-size', '10px').text(`${d.completeness}%`);
    centerChartGroup.append('text').attr('x', -barWidth / 2 - barGap / 2).attr('y', maxBarHeight + 12).attr('text-anchor', 'middle').attr('font-size', '10px').attr('font-weight', 'bold').text('Comp.');
    const contScaleMax = actualMaxContDisplay > 0 ? actualMaxContDisplay : 1;
    const contHeight = (d.contamination / contScaleMax) * maxBarHeight;
    centerChartGroup.append('rect').attr('x', barGap / 2).attr('y', maxBarHeight - contHeight).attr('width', barWidth).attr('height', contHeight).attr('class', 'modal-contamination-bar').attr('rx', 3);
    centerChartGroup.append('text').attr('x', barWidth / 2 + barGap / 2).attr('y', maxBarHeight - contHeight - 5).attr('text-anchor', 'middle').attr('font-size', '10px').text(`${d.contamination}%`);
    centerChartGroup.append('text').attr('x', barWidth / 2 + barGap / 2).attr('y', maxBarHeight + 12).attr('text-anchor', 'middle').attr('font-size', '10px').attr('font-weight', 'bold').text('Cont.');
  }
}

// Smart tooltip positioning function
function positionTooltip(event, tooltipContent) {
  tooltip.html(tooltipContent);

  // Get tooltip dimensions
  const tooltipNode = tooltip.node();
  const tooltipRect = tooltipNode.getBoundingClientRect();
  const tooltipWidth = tooltipRect.width;
  const tooltipHeight = tooltipRect.height;

  // Get viewport dimensions
  const viewportWidth = window.innerWidth;
  const viewportHeight = window.innerHeight;

  // Default offset from cursor
  const defaultOffsetX = 10;
  const defaultOffsetY = 10;

  // Calculate potential positions
  let left = event.pageX + defaultOffsetX;
  let top = event.pageY + defaultOffsetY;

  // Check if tooltip would go off the right edge
  if (left + tooltipWidth > viewportWidth) {
    left = event.pageX - tooltipWidth - defaultOffsetX; // Position to the left of cursor
  }

  // Check if tooltip would go off the bottom edge
  if (top + tooltipHeight > viewportHeight) {
    top = event.pageY - tooltipHeight - defaultOffsetY; // Position above cursor
  }

  // Ensure tooltip doesn't go off the left edge
  if (left < 0) {
    left = defaultOffsetX;
  }

  // Ensure tooltip doesn't go off the top edge
  if (top < 0) {
    top = defaultOffsetY;
  }

  tooltip.style('left', `${left}px`).style('top', `${top}px`).style('opacity', 1);
}

// Initial setup and data loading
d3.json("visualization_data.json?v=" + Date.now()).then(function(loadedData) {
  vizData = loadedData;

  // Populate sample dropdown
  sampleDropdown
    .selectAll('option')
    .data(Object.keys(vizData))
    .enter()
    .append('option')
    .attr('value', d => d)
    .text(d => d);
  sampleDropdown.property('value', Object.keys(vizData)[0] || '');

  // Setup legend item click handlers
  d3.selectAll('.legend-item').on('click', function() {
    const quality = d3.select(this).attr('data-quality');
    if (activeQualityFilters.has(quality)) {
      activeQualityFilters.delete(quality);
      d3.select(this).classed('filter-active', false);
    } else {
      activeQualityFilters.add(quality);
      d3.select(this).classed('filter-active', true);
    }
    renderChart(); // Re-render the chart with new filters
  });

  // Setup new view toggle logic
  const viewOptions = d3.selectAll('#details-view-toggle-container .view-option');
  viewOptions.on('click', function() {
    viewOptions.classed('active', false);
    d3.select(this).classed('active', true);

    // Store the selected view type (e.g., 'circular', 'bar', 'gc')
    // This will be read by updateDetailsView
    window.selectedDetailsView = d3.select(this).attr('data-view');

    const currentData = window.currentMagData;
    const currentMaxCont = window.currentMaxContDisplay;
    if (currentData) {
      updateDetailsView(currentData, currentMaxCont);
    }
  });

  // Set initial view state
  window.selectedDetailsView = 'circular'; // Default to circular

  // Remove the old toggle switch if it exists (no longer needed)
  const oldViewToggleElement = d3.select('#view-toggle'); // Select the old element by ID
  if (!oldViewToggleElement.empty()){
     oldViewToggleElement.remove(); // Remove it if it exists
  }

  // Setup taxonomy distribution modal
  const taxonomyBtn = d3.select('#taxonomy-distribution-btn');
  const taxonomyOverlay = d3.select('#taxonomy-overlay');
  const taxonomyModalClose = d3.select('#taxonomy-modal-close');

  taxonomyBtn.on('click', function() {
    showTaxonomyDistribution();
  });

  taxonomyModalClose.on('click', function() {
    taxonomyOverlay.style('display', 'none');
  });

  taxonomyOverlay.on('click', function(event) {
    if (event.target === taxonomyOverlay.node()) {
      taxonomyOverlay.style('display', 'none');
    }
  });

  renderChart(); // Initial render
  updateContaminationValueDisplay(); // Initial call to set display value
  updateCompletenessValueDisplay(); // Initial call for completeness display
  d3.select("#collapse-level-select").property("value", "none"); // Ensure default is no collapse

  // Setup reset filters button
  const resetFiltersBtn = d3.select('#reset-filters-btn');
  if (!resetFiltersBtn.empty()) {
    resetFiltersBtn.on('click', function() {
      // Reset contamination slider
      contThresholdSlider.property('value', 100);
      updateContaminationValueDisplay();

      // Reset completeness slider
      compThresholdSlider.property('value', 0);
      updateCompletenessValueDisplay();

      // Reset quality filters
      activeQualityFilters.clear();
      d3.selectAll('.legend-item').classed('filter-active', false);

      // Reset other dropdowns to their default visual state and value
      d3.select('#sort-select').property('value', 'completeness');
      d3.select('#order-select').property('value', 'descending');
      d3.select('#scaling-select').property('value', 'length');
      d3.select('#classification-level-select').property('value', 'none');
      d3.select('#collapse-level-select').property('value', 'none');
      
      // Note: Sample dropdown is not reset by this button, 
      // as filter reset usually applies to the current sample view.

      renderChart(); // Re-render the chart with reset filters
    });
  }

}).catch(function(error) {
  console.error("Error loading the data: ", error);
  warning.text("Error loading visualization data. Please check the console.").style('display', 'block');
});

// Function to show taxonomy distribution modal
function showTaxonomyDistribution() {
  const selectedSample = sampleDropdown.property('value');
  const rawMags = vizData[selectedSample] || [];
  // Suggested log: console.log('Raw MAGs for sample '+ selectedSample + ':', rawMags.length, rawMags);

  if (rawMags.length === 0) {
    alert('No MAGs found in the selected sample.');
    return;
  }
  
  let filteredMags = [...rawMags]; // Start with a copy of rawMags

  // Apply quality filters
  if (activeQualityFilters.size > 0) {
    filteredMags = filteredMags.filter(mag => activeQualityFilters.has(mag.quality));
  }
  // Suggested log: console.log('After quality filters ('+ activeQualityFilters.size + ' active):', filteredMags.length, activeQualityFilters);

  // Apply classification level filter
  const selectedClassificationLevel = d3.select('#classification-level-select').property('value');
  filteredMags = filteredMags.filter(mag => 
    meetsClassificationLevel(mag.assigned_taxonomy, selectedClassificationLevel)
  );
  // Suggested log: console.log('After classification filter ('+ selectedClassificationLevel + '):', filteredMags.length);
  
  // Apply contamination filter
  const contaminationFilter = +d3.select('#cont-threshold').property('value');
  filteredMags = filteredMags.filter(d => d.contamination <= contaminationFilter);
  // Suggested log: console.log('After contamination filter (<='+ contaminationFilter + '%):', filteredMags.length);
  
  // Apply completeness filter
  const completenessFilter = +d3.select('#comp-threshold').property('value');
  filteredMags = filteredMags.filter(d => d.completeness >= completenessFilter);
  // Suggested log: console.log('After completeness filter (>='+ completenessFilter + '%):', filteredMags.length);

  if (filteredMags.length === 0) {
    alert('No MAGs match the current filters.');
    return;
  }
  // Suggested log: console.log('Final filtered MAGs for taxonomy chart:', filteredMags.length, filteredMags);
  
  d3.select('#taxonomy-modal-title').text(`Taxonomy Distribution - ${selectedSample} (${filteredMags.length} MAGs)`);
  d3.select('#taxonomy-overlay').style('display', 'flex');
  renderTaxonomyDistributionChart(filteredMags);
}

// Function to render taxonomy distribution chart
function renderTaxonomyDistributionChart(mags) {
  // Suggested log: console.log('renderTaxonomyDistributionChart received MAGs:', mags.length, mags);
  const svg = d3.select('#taxonomy-chart');
  svg.selectAll('*').remove();

  // Count MAGs by taxonomic level
  const levelCounts = {
    'Level 1': 0,
    'Level 2': 0,
    'Level 3': 0,
    'Level 4': 0,
    'Level 5': 0,
    'Level 6': 0,
    'Level 7': 0,
    'Level 8': 0,
    'Unknown': 0
  };

  mags.forEach(mag => {
    const deepestLevel = getDeepestTaxonomicLevel(mag.assigned_taxonomy);
    if (!deepestLevel) {
      levelCounts['Unknown']++;
    } else {
      const levelMap = {
        'd': 'Level 1',
        'k': 'Level 2',
        'p': 'Level 3',
        'c': 'Level 4',
        'o': 'Level 5',
        'f': 'Level 6',
        'g': 'Level 7',
        's': 'Level 8'
      };
      const levelName = levelMap[deepestLevel] || 'Unknown';
      levelCounts[levelName]++;
    }
  });
  // Suggested log: console.log('Taxonomy levelCounts:', levelCounts);
  
  // Convert to array and filter out zeros
  const data = Object.entries(levelCounts)
    .filter(([level, count]) => count > 0)
    .map(([level, count]) => ({ level, count }))
    .sort((a, b) => {
      // Sort by level number, with Unknown at the end
      if (a.level === 'Unknown') return 1;
      if (b.level === 'Unknown') return -1;
      const aNum = parseInt(a.level.split(' ')[1]);
      const bNum = parseInt(b.level.split(' ')[1]);
      return aNum - bNum;
    });
  // Suggested log: console.log('Data for taxonomy bars:', data);
  
  if (data.length === 0) {
    svg.append('text')
      .attr('x', '50%')
      .attr('y', '50%')
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'central')
      .text('No taxonomy data found');
    return;
  }

  // Set up dimensions
  const containerNode = svg.node().parentNode;
  const containerWidth = containerNode.clientWidth || 550;
  const containerHeight = containerNode.clientHeight || 300;
  const margin = { top: 20, right: 20, bottom: 60, left: 80 };
  const width = containerWidth - margin.left - margin.right;
  const height = containerHeight - margin.top - margin.bottom;

  svg.attr('width', containerWidth)
     .attr('height', containerHeight);

  // Create scales
  const xScale = d3.scaleBand()
    .domain(data.map(d => d.level))
    .range([0, width])
    .padding(0.2);

  const yScale = d3.scaleLinear()
    .domain([0, d3.max(data, d => d.count)])
    .range([height, 0]);

  const colorScale = d3.scaleSequential()
    .domain([0, data.length - 1])
    .interpolator(d3.interpolateBlues);

  // Create main group
  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  // Add bars
  g.selectAll('.bar')
    .data(data)
    .enter()
    .append('rect')
    .attr('class', 'bar')
    .attr('x', d => xScale(d.level))
    .attr('y', d => yScale(d.count))
    .attr('width', xScale.bandwidth())
    .attr('height', d => height - yScale(d.count))
    .attr('fill', (d, i) => colorScale(i))
    .attr('rx', 4)
    .on('mouseover', function(event, d) {
      d3.select(this).style('opacity', 0.8);
      tooltip.style('opacity', 1);
    })
    .on('mousemove', function(event, d) {
      const percentage = ((d.count / mags.length) * 100).toFixed(1);
      positionTooltip(event, `<strong>${d.level}</strong><br/>Count: ${d.count}<br/>Percentage: ${percentage}%`);
    })
    .on('mouseout', function() {
      d3.select(this).style('opacity', 1);
      tooltip.style('opacity', 0);
    });

  // Add value labels on top of bars
  g.selectAll('.bar-label')
    .data(data)
    .enter()
    .append('text')
    .attr('class', 'bar-label')
    .attr('x', d => xScale(d.level) + xScale.bandwidth() / 2)
    .attr('y', d => yScale(d.count) - 5)
    .attr('text-anchor', 'middle')
    .attr('font-size', '12px')
    .attr('font-weight', '600')
    .attr('fill', 'var(--text-primary)')
    .text(d => d.count);

  // Add x-axis
  g.append('g')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(xScale))
    .selectAll('text')
    .style('text-anchor', 'middle')
    .attr('font-size', '11px');

  // Add y-axis
  g.append('g')
    .call(d3.axisLeft(yScale).ticks(Math.min(10, d3.max(data, d => d.count))))
    .selectAll('text')
    .attr('font-size', '11px');

  // Add axis labels
  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('y', 0 - margin.left)
    .attr('x', 0 - (height / 2))
    .attr('dy', '1em')
    .style('text-anchor', 'middle')
    .attr('font-size', '12px')
    .attr('font-weight', '600')
    .text('Number of MAGs');

  g.append('text')
    .attr('transform', `translate(${width / 2}, ${height + margin.bottom - 10})`)
    .style('text-anchor', 'middle')
    .attr('font-size', '12px')
    .attr('font-weight', '600')
    .text('Deepest Classification Level');
}

// Update event listeners
d3.selectAll('#sample-select,#sort-select,#order-select,#scaling-select,#classification-level-select,#collapse-level-select').on('change', renderChart);

contThresholdSlider.on('input', function() {
renderChart();
  updateContaminationValueDisplay();
});

compThresholdSlider.on('input', function() {
  renderChart();
  updateCompletenessValueDisplay();
});

// Function to render bar chart view in modal
function renderBarChartView(d, actualMaxContDisplay) {
  const svgNode = detailsSVG.node();
  const w = svgNode.clientWidth || 300; // Use actual width with fallback
  const h = svgNode.clientHeight || 320; // Use actual height with fallback
  const margin = { top: 25, right: 12, bottom: 10, left: 12 }; // Increased margins further
  const chartWidth = w - margin.left - margin.right;
  const availableHeight = h - margin.top - margin.bottom;

  const contigsMap = d.contigs || {};

  // Calculate total count per taxon and sort by count (descending)
  const taxonData = Object.entries(contigsMap).map(([groupTaxon, specificContigsArray]) => ({
    taxon: groupTaxon,
    count: specificContigsArray.length,
    totalLength: d3.sum(specificContigsArray, c => c.length)
  })).sort((a, b) => b.count - a.count);

  if (taxonData.length === 0) {
    detailsSVG.append("text")
      .attr("x", w / 2).attr("y", h / 2)
      .attr("text-anchor", "middle").attr("dominant-baseline", "central")
      .text("No contigs to display.");
    return;
  }

  // Create color scale based on taxa ordered by count (same as main chart)
  const orderedTaxa = taxonData.map(d => d.taxon);
  const colorScale = makeColorScale(orderedTaxa);

  // Fixed dimensions for consistent appearance - maximize space usage
  const rowHeight = 20;
  const rowGap = 6;
  const maxBarWidth = chartWidth * 0.65; // Limited to 50% for comfortable label spacing

  // Calculate total height needed for all bars
  const totalContentHeight = taxonData.length * (rowHeight + rowGap);
  const needsScroll = totalContentHeight > availableHeight;

  // Set SVG viewBox to accommodate all content, but don't change the actual SVG size
  const svgContentHeight = totalContentHeight + margin.top + margin.bottom;
  detailsSVG.attr('viewBox', `0 0 ${w} ${svgContentHeight}`)
            .attr('preserveAspectRatio', 'xMidYMin meet');

  // Don't change the parent container styling - let CSS handle the scrolling
  // The overflow: auto on #details-chart will handle scrolling automatically

  // Create scales
  const yScale = d3.scaleBand()
    .domain(orderedTaxa)
    .range([margin.top, margin.top + totalContentHeight])
    .padding(0);

  const xScale = d3.scaleLinear()
    .domain([0, d3.max(taxonData, d => d.count)])
    .range([0, maxBarWidth]);

  const chartGroup = detailsSVG.append('g');

  // Draw horizontal bars starting from left edge
  chartGroup.selectAll('.horizontal-bar')
    .data(taxonData)
    .enter()
    .append('rect')
    .attr('class', 'horizontal-bar')
    .attr('x', margin.left)
    .attr('y', d => yScale(d.taxon))
    .attr('width', d => xScale(d.count))
    .attr('height', rowHeight)
    .attr('fill', d => colorScale(d.taxon))
    .attr('rx', 5)
    .on('mouseover', function(event, d) {
      d3.select(this).style('fill', d3.color(colorScale(d.taxon)).brighter(0.7));
      tooltip.style('opacity', 1);
    })
    .on('mousemove', function(event, d) {
      const avgLength = Math.round(d.totalLength / d.count);
      const tooltipContent = `<strong>Taxon Group</strong><br/>` +
                            `Group: ${parseUserFriendlyTaxonomy(d.taxon)}<br/>` +
                            `Contig Count: ${d.count}<br/>` +
                            `Total Length: ${d.totalLength.toLocaleString()} bp<br/>` +
                            `Average Length: ${avgLength.toLocaleString()} bp` +
                            `${formatDetailedTaxonomy(d.taxon)}`;
      positionTooltip(event, tooltipContent);
    })
    .on('mouseout', function(event, d) {
      d3.select(this).style('fill', colorScale(d.taxon));
      tooltip.style('opacity', 0);
    });

  // Add count labels inside bars at the left edge (only if bar is wide enough)
  chartGroup.selectAll('.count-label')
    .data(taxonData.filter(d => xScale(d.count) > 18)) // Adjusted threshold
    .enter()
    .append('text')
    .attr('class', 'count-label')
    .attr('x', d => margin.left + 6) // Position at left edge of bars with small padding
    .attr('y', d => yScale(d.taxon) + rowHeight / 2)
    .attr('dy', '0.35em')
    .attr('text-anchor', 'start') // Changed from 'end' to 'start'
    .attr('font-size', '9px')
    .attr('font-weight', '600')
    .attr('fill', 'white')
    .style('text-shadow', '1px 1px 1px rgba(0,0,0,0.5)')
    .text(d => d.count);

  // Add taxon labels on the right side
  chartGroup.selectAll('.taxon-label')
    .data(taxonData)
    .enter()
    .append('text')
    .attr('class', 'taxon-label')
    .attr('x', d => margin.left + xScale(d.count) + 8) // Position right after bar end with small gap
    .attr('y', d => yScale(d.taxon) + rowHeight / 2)
    .attr('dy', '0.35em')
    .attr('font-size', '10px')
    .attr('fill', 'var(--text-primary)')
    .text(d => {
      const friendlyName = parseUserFriendlyTaxonomy(d.taxon);
      // Calculate remaining space from bar end to panel edge
      const remainingSpace = w - margin.right - (margin.left + xScale(d.count) + 8);
      const maxChars = Math.floor(remainingSpace / 5.5);
      return friendlyName.length > maxChars ?
             friendlyName.substring(0, maxChars - 2) + '...' :
             friendlyName;
    });

  // Add count labels next to narrow bars on the left side
  chartGroup.selectAll('.external-count-label')
    .data(taxonData.filter(d => xScale(d.count) <= 18)) // Show external count for narrow bars
    .enter()
    .append('text')
    .attr('class', 'external-count-label')
    .attr('x', margin.left - 3) // Position just before the bar starts
    .attr('y', d => yScale(d.taxon) + rowHeight / 2)
    .attr('dy', '0.35em')
    .attr('text-anchor', 'end') // Right-align so it doesn't overlap with bars
    .attr('font-size', '9px')
    .attr('font-weight', '600')
    .attr('fill', 'var(--text-secondary)')
    .text(d => d.count); // Remove parentheses since it's clearly separated now

  // Add title
  chartGroup.append('text')
    .attr('x', w / 2)
    .attr('y', margin.top - 8) // Closer to top
    .attr('text-anchor', 'middle')
    .attr('font-size', '12px')
    .attr('font-weight', 'bold')
    .attr('fill', 'var(--text-primary)')
    .text('Contigs per Taxon');

  // Add scroll indicator if needed
  if (needsScroll) {
    chartGroup.append('text')
      .attr('x', w - margin.right - 1)
      .attr('y', margin.top - 2)
      .attr('text-anchor', 'end')
      .attr('font-size', '8px')
      .attr('fill', 'var(--text-secondary)')
      .style('opacity', 0.7)
      .text('↕ scroll');
  }
}

// Helper function to calculate length-weighted or unweighted quantiles
function calculateQuantile(contigsData, quantile, isWeighted) {
  if (!contigsData || contigsData.length === 0) {
    return undefined;
  }

  // Ensure contigs are sorted by gc_content
  const sortedContigs = [...contigsData].sort((a, b) => a.gc_content - b.gc_content);

  const totalWeight = isWeighted ? d3.sum(sortedContigs, c => c.length) : sortedContigs.length;
  if (totalWeight === 0) return undefined;

  const targetWeight = totalWeight * quantile;
  let cumulativeWeight = 0;

  for (const contig of sortedContigs) {
    const currentWeight = isWeighted ? contig.length : 1;
    if (cumulativeWeight + currentWeight >= targetWeight) {
      return contig.gc_content;
    }
    cumulativeWeight += currentWeight;
  }
  return sortedContigs[sortedContigs.length - 1].gc_content; // Fallback
}

// Function to render GC content distribution as KDE plot in modal
function renderMAGGCKDEPlot(d) {
  detailsSVG.selectAll('*').remove(); // Clear SVG immediately

  window.requestAnimationFrame(() => {
    // Calculate stable display dimensions for the SVG element
    const detailsColumnNode = d3.select('#mag-details-column').node();
    const titleNode = d3.select('#details-title').node();
    const metricsNode = d3.select('#details-metrics-text').node();

    const titleActualHeight = titleNode.offsetHeight + parseFloat(getComputedStyle(titleNode).marginBottom || '0');
    const metricsActualHeight = metricsNode.offsetHeight + parseFloat(getComputedStyle(metricsNode).marginBottom || '0');

    let svgDisplayHeight = (detailsColumnNode.clientHeight || 320) - titleActualHeight - metricsActualHeight;
    const cssMinHeightSVG = 280; // From #details-chart min-height in CSS
    if (svgDisplayHeight < cssMinHeightSVG) {
      svgDisplayHeight = cssMinHeightSVG;
    }
    const svgDisplayWidth = detailsColumnNode.clientWidth > 0 ? detailsColumnNode.clientWidth : 300;

    // Explicitly set SVG element's actual style width and height
    detailsSVG
      .style('width', svgDisplayWidth + 'px')
      .style('height', svgDisplayHeight + 'px');

    const margin = { top: 30, right: 20, bottom: 50, left: 70 };
    const fixedPlotAreaHeight = 200; // Define a fixed height for the main plot drawing area
    const checkboxHeightWithPadding = 30;

    // Width for the chart's internal drawing area (within margins)
    const chartDrawingWidth = svgDisplayWidth - margin.left - margin.right;

    if (chartDrawingWidth <= 0) {
      detailsSVG.append("text").attr("x", svgDisplayWidth / 2).attr("y", svgDisplayHeight / 2).attr("text-anchor", "middle").text("Container too small for chart.");
      return;
    }

    // Calculate total internal content height for the viewBox
    const totalInternalViewBoxHeight = margin.top + fixedPlotAreaHeight + margin.bottom + checkboxHeightWithPadding;

    // Set viewBox using the SVG's display width and the calculated total internal content height
    detailsSVG
      .attr('viewBox', `0 0 ${svgDisplayWidth} ${totalInternalViewBoxHeight}`)
      .attr('preserveAspectRatio', 'xMidYMin meet');

    const contigsMap = d.contigs || {};
    const allMAGContigs = [];
    for (const [_groupTaxon, specificContigsArray] of Object.entries(contigsMap)) {
      if (!specificContigsArray || specificContigsArray.length === 0) continue;
      specificContigsArray.forEach(c => {
        if (c.gc_content !== undefined && c.gc_content !== null && c.length !== undefined && c.length > 0) {
          allMAGContigs.push({ gc_content: c.gc_content, length: c.length });
        }
      });
    }

    if (allMAGContigs.length === 0) {
      detailsSVG.append("text")
        .attr("x", svgDisplayWidth / 2) // Centered based on SVG display width
        .attr("y", margin.top + fixedPlotAreaHeight / 2) // Centered in the fixed plot area
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "central")
        .text("No GC content data to display for this MAG.");
      return;
    }

    // Calculate Length-Weighted KDE
    function calculateKDE(data, bandwidth, isWeightedCalc) {
      const kdePoints = [];
      const xDomain = d3.range(0, 100.5, 0.5); // GC points to evaluate
      let maxDensity = 0;

      for (const x_eval of xDomain) {
        let density_at_x_eval = 0;
        for (const contig of data) {
          const u = (x_eval - contig.gc_content) / bandwidth;
          const kernel_value = (1 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * u * u); // Gaussian kernel
          const weight = isWeightedCalc ? contig.length : 1;
          density_at_x_eval += kernel_value * weight;
        }
        kdePoints.push({gc: x_eval, density: density_at_x_eval });
        if (density_at_x_eval > maxDensity) maxDensity = density_at_x_eval;
      }
      return { points: kdePoints, maxDensity: maxDensity };
    }

    const bandwidth = 2.5; // KDE bandwidth - adjust for smoothness
    const { points: kdeData, maxDensity: calculatedMaxDensity } = calculateKDE(allMAGContigs, bandwidth, window.gcKDEIsWeighted);

    // Calculate length-weighted or unweighted mean and median for the entire MAG
    let totalWeightSum = 0;
    let weightedGCSum = 0;
    allMAGContigs.forEach(c => {
      const currentWeight = window.gcKDEIsWeighted ? c.length : 1;
      totalWeightSum += currentWeight;
      weightedGCSum += c.gc_content * currentWeight;
    });
    const meanGC = totalWeightSum > 0 ? weightedGCSum / totalWeightSum : undefined;
    const medianGC = calculateQuantile(allMAGContigs, 0.50, window.gcKDEIsWeighted);

    // Create Scales
    const xScale = d3.scaleLinear()
      .domain([0, 100]) // GC content from 0 to 100%
      .range([0, chartDrawingWidth]); // USE chartDrawingWidth

    const yScale = d3.scaleLinear()
      .domain([0, calculatedMaxDensity > 0 ? calculatedMaxDensity * 1.1 : 1]) // Add 10% headroom or default
      .range([fixedPlotAreaHeight, 0]); // USE fixedPlotAreaHeight

    const chartGroup = detailsSVG.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Draw KDE area and line
    const area = d3.area()
      .x(d_kde => xScale(d_kde.gc))
      .y0(fixedPlotAreaHeight) // USE fixedPlotAreaHeight
      .y1(d_kde => yScale(d_kde.density))
      .curve(d3.curveBasis); // Smooth curve

    chartGroup.append('path')
      .datum(kdeData)
      .attr('class', 'kde-area')
      .attr('d', area)
      .style('fill', 'var(--accent-primary)')
      .style('opacity', 0.3);

    const line = d3.line()
      .x(d_kde => xScale(d_kde.gc))
      .y(d_kde => yScale(d_kde.density))
      .curve(d3.curveBasis);

    chartGroup.append('path')
      .datum(kdeData)
      .attr('class', 'kde-line')
      .attr('d', line)
      .style('stroke', 'var(--accent-primary)')
      .style('stroke-width', 2)
      .style('fill', 'none');

    // Draw vertical lines for mean and median
    if (meanGC !== undefined) {
      chartGroup.append('line')
        .attr('x1', xScale(meanGC))
        .attr('x2', xScale(meanGC))
        .attr('y1', fixedPlotAreaHeight) // USE fixedPlotAreaHeight
        .attr('y2', 0)
        .attr('stroke', 'var(--accent-secondary)')
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '4,2');
      chartGroup.append('text')
          .attr('x', xScale(meanGC) + 4)
          .attr('y', 10)
          .attr('fill', 'var(--accent-secondary)')
          .style('font-size', '9px')
          .text(`Mean: ${meanGC.toFixed(1)}%`);
    }
    if (medianGC !== undefined) {
      chartGroup.append('line')
        .attr('x1', xScale(medianGC))
        .attr('x2', xScale(medianGC))
        .attr('y1', fixedPlotAreaHeight) // USE fixedPlotAreaHeight
        .attr('y2', 0)
        .attr('stroke', 'var(--danger-color)') // Different color for median
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '4,2');
       chartGroup.append('text')
          .attr('x', xScale(medianGC) + 4)
          .attr('y', 22)
          .attr('fill', 'var(--danger-color)')
          .style('font-size', '9px')
          .text(`Median: ${medianGC.toFixed(1)}%`);
    }

    // Add Axes
    chartGroup.append('g').attr('class', 'axis').call(d3.axisLeft(yScale).ticks(5).tickFormat(d3.format(".2s"))); // Format Y for density values
    chartGroup.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('y', 0 - margin.left + 15)
        .attr('x', 0 - (fixedPlotAreaHeight / 2)) // USE fixedPlotAreaHeight
        .attr('dy', '0em')
        .attr('class', 'chart-axis-label')
        .style('text-anchor', 'middle')
        .text(window.gcKDEIsWeighted ? 'Length-Weighted Density' : 'Density');

    chartGroup.append('g')
      .attr('transform', `translate(0,${fixedPlotAreaHeight})`) // USE fixedPlotAreaHeight
      .attr('class', 'axis')
      .call(d3.axisBottom(xScale).ticks(10));
    chartGroup.append('text')
        .attr('x', chartDrawingWidth / 2) // USE chartDrawingWidth
        .attr('y', fixedPlotAreaHeight + margin.bottom - 10) // USE fixedPlotAreaHeight
        .attr('class', 'chart-axis-label')
        .style('text-anchor', 'middle')
        .text('GC Content (%)');

    // Tooltip for KDE plot
    const focus = chartGroup.append("g")
        .attr("class", "focus")
        .style("display", "none");

    focus.append("line") // Vertical line following mouse
        .attr("class", "focus-line")
        .attr("y1", 0)
        .attr("y2", fixedPlotAreaHeight) // USE fixedPlotAreaHeight
        .style("stroke", "#333")
        .style("stroke-width", 1)
        .style("stroke-dasharray", "3,3");

    detailsSVG.append('rect') // Overlay for mouse events
      .attr('transform', `translate(${margin.left},${margin.top})`)
      .attr('class', 'overlay-rect')
      .attr('width', chartDrawingWidth) // USE chartDrawingWidth
      .attr('height', fixedPlotAreaHeight) // USE fixedPlotAreaHeight
      .style('fill', 'none')
      .style('pointer-events', 'all')
      .on('mouseover', () => { focus.style("display", null); tooltip.style('opacity', 1); })
      .on('mouseout', () => { focus.style("display", "none"); tooltip.style('opacity', 0); })
      .on('mousemove', (event) => {
        const x_mouse_relative_to_overlay = d3.pointer(event, event.currentTarget)[0];
        const x0 = xScale.invert(x_mouse_relative_to_overlay);

        const bisectGC = d3.bisector(d_kde => d_kde.gc).left; // Renamed for clarity
        const i = bisectGC(kdeData, x0, 1);

        const d0 = kdeData[i - 1];
        const d1 = kdeData[i];

        let d_point;
        if (d0 && d1) { // Standard case: mouse is between two points
          d_point = (x0 - d0.gc > d1.gc - x0) ? d1 : d0;
        } else if (d0) { // Mouse is at or beyond the last data point
          d_point = d0;
        } else if (d1) { // Mouse is at or before the first data point
          d_point = d1;
        } else { // Should not happen if kdeData is not empty
          return;
        }

        if (!d_point) { // Fallback, if something unexpected occurred
            return;
        }

        focus.select(".focus-line").attr("transform", `translate(${xScale(d_point.gc)},0)`);
        const meanMedianLabel = window.gcKDEIsWeighted ? "LW" : "Unweighted";
        positionTooltip(event,
          `GC: ${d_point.gc.toFixed(1)}%<br/>` +
          `Density: ${d_point.density.toExponential(2)}<br/>` +
          (meanGC ? `${meanMedianLabel} Mean: ${meanGC.toFixed(1)}%<br/>` : '') +
          (medianGC ? `${meanMedianLabel} Median: ${medianGC.toFixed(1)}%` : '')
        );
      });

    // Add chart title and toggle
    const titleText = "MAG GC Distribution"; // Removed parentheses
    detailsSVG.append("text")
      .attr("x", svgDisplayWidth / 2)
      .attr("y", margin.top / 2 + 2) // Adjusted y to make space for toggle if needed
      .attr("text-anchor", "middle")
      .style("font-size", "12px")
      .style("font-weight", "bold")
      .text(titleText);

    // Add Length-Weighted Toggle Checkbox
    const foreignObject = detailsSVG.append("foreignObject")
      .attr("x", svgDisplayWidth / 2 - 65)
      .attr("y", margin.top + fixedPlotAreaHeight + margin.bottom + 10) // USE fixedPlotAreaHeight
      .attr("width", 130)
      .attr("height", 20);

    const checkboxLabel = foreignObject.append("xhtml:label")
      .style("font-size", "10px")
      .style("font-family", "var(--font-sans-serif)")
      .style("color", "var(--text-secondary)")
      .style("display", "flex")
      .style("align-items", "center");

    checkboxLabel.append("xhtml:input")
      .attr("type", "checkbox")
      .property("checked", window.gcKDEIsWeighted)
      .style("margin-right", "4px")
      .on("change", function() {
        window.gcKDEIsWeighted = this.checked;
        renderMAGGCKDEPlot(d); // Re-render with new setting
      });
    checkboxLabel.append("xhtml:span")
      .text("Length-Weighted");
  }); // End of requestAnimationFrame callback
}

// --- Info Icon Tooltip Logic ---
const mainChartInfoIcon = d3.select('#main-chart-info-icon');
const detailsPanelInfoIcon = d3.select('#details-panel-info-icon');

const mainChartExplanationText = `
  <strong>Main MAGs Chart:</strong><br/>
  This chart visualizes Metagenome-Assembled Genomes (MAGs) for the selected sample. Each row is a MAG.<br/><br/>
  <strong>Left Side (Quality Metrics):</strong><br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;The dark grey bar (left of center) shows MAG <strong>Contamination (%)</strong>. Its length is relative to the maximum contamination in the current filtered set (or 1% if max is 0).<br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;The green bar (right of center) shows MAG <strong>Completeness (%)</strong> on a 0-100% scale.<br/><br/>
  <strong>Right Side (Assembly Composition):</strong><br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;This stacked bar shows MAG composition. Scale by <strong>Assembly Length</strong> (total bp) or <strong>Contig Count</strong> via dropdown.<br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;Segments are taxonomic groups. Colors are based on dominant taxa (by length) across all MAGs. Hover for details.<br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;Use 'Collapse by' to aggregate to a higher taxonomic level. Color is from the dominant original taxon in the collapsed segment.<br/><br/>
  <strong>Far Right (Quality Tier):</strong><br/>
  &nbsp;&nbsp;&bull;&nbsp;&nbsp;Icon indicates MAG quality (High <i class='bi bi-check-circle-fill' style='color:var(--success-color); vertical-align: -0.1em;'></i>, Medium <i class='bi bi-circle-half' style='color:var(--accent-primary); vertical-align: -0.1em;'></i>, Low <i class='bi bi-exclamation-triangle-fill' style='color:var(--danger-color); vertical-align: -0.1em;'></i>).<br/><br/>
  <em>Click a MAG row to see details in the right panel.</em>
`;

function getDetailsPanelExplanation() {
  const magId = window.currentMagData ? window.currentMagData.id : "Selected MAG";
  let explanation = `<strong>MAG Details Panel (${magId}):</strong><br/>`;
  explanation += "Switch views using 'Circular', 'Bar Chart', 'GC Dist.' in controls.<br/><br/>";

  const currentView = window.selectedDetailsView || 'circular';

  if (currentView === 'circular') {
    explanation += `
      <strong>Circular Contig Plot:</strong><br/>
      Displays contigs of the MAG. The full circle represents 100% completeness; plotted arcs span a portion proportional to actual completeness.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;<strong>Inner Ring (Contigs):</strong> Colored arcs are contigs or aggregated small contigs from the same taxon group. Arc angle is proportional to its length relative to total MAG contig length (within completeness-scaled circle). Colors match main chart.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;<strong>Outer Ring (Coverage):</strong> Grey segments correspond to inner contigs. Radial thickness is proportional to average contig coverage, scaled to max coverage in this MAG.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;<strong>Center:</strong> Bar chart shows overall MAG Completeness and Contamination.<br/>
      <em>Hover over arcs for details.</em>`;
  } else if (currentView === 'bar') {
    explanation += `
      <strong>Contigs per Taxon (Bar Chart):</strong><br/>
      Horizontal bars show the number of contigs for each taxon group in this MAG.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;Bars are ordered by contig count (descending). Colors match main chart.<br/>
      <em>Hover over bars for details.</em>`;
  } else if (currentView === 'gc') {
    explanation += `
      <strong>GC Content Distribution (KDE):</strong><br/>
      A Kernel Density Estimate plot of GC content for all contigs in this MAG.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;X-axis: GC content (%). Y-axis: Density.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;Vertical lines: Mean and Median GC for the MAG.<br/>
      &nbsp;&nbsp;&bull;&nbsp;&nbsp;Use 'Length-Weighted' checkbox below to toggle between unweighted KDE (contigs count equally) and length-weighted KDE (contigs weighted by length).<br/>
      <em>Hover over curve for GC/density.</em>`;
  } else {
    explanation += "Select a view to see details.";
  }
  return explanation;
}

if (!mainChartInfoIcon.empty()) {
  mainChartInfoIcon
    .on('mouseover', () => tooltip.style('opacity', 1))
    .on('mousemove', (event) => positionTooltip(event, mainChartExplanationText))
    .on('mouseout', () => tooltip.style('opacity', 0));
}

if (!detailsPanelInfoIcon.empty()) {
  detailsPanelInfoIcon
    .on('mouseover', () => tooltip.style('opacity', 1))
    .on('mousemove', (event) => positionTooltip(event, getDetailsPanelExplanation()))
    .on('mouseout', () => tooltip.style('opacity', 0));
}
// --- End Info Icon Tooltip Logic ---
