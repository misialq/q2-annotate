"use strict";

const tfa = JSON.parse(document.getElementById("tfa-data").textContent);
const taxonSelect = document.getElementById("taxon");
const functionSelect = document.getElementById("function");
const groupSelect = document.getElementById("group");
const limitSelect = document.getElementById("limit");
const status = document.getElementById("status");
const schema = "https://vega.github.io/schema/vega-lite/v5.json";
let activeViews = [];
let renderQueue = Promise.resolve();

function addOptions(select, labels, firstLabel) {
  const first = document.createElement("option");
  first.value = "";
  first.textContent = firstLabel;
  select.appendChild(first);
  for (const label of labels) {
    const option = document.createElement("option");
    option.value = label;
    option.textContent = label;
    select.appendChild(option);
  }
}

const taxa = [...new Set(tfa.pairs.map(pair => pair.taxon))].sort();
const functions = [...new Set(tfa.pairs.map(pair => pair.function))].sort();
const loadsByPair = tfa.pairs.map(() => []);
for (const load of tfa.loads) loadsByPair[load[0]].push(load);
addOptions(taxonSelect, taxa, "All taxa");
addOptions(functionSelect, functions, "All functions");
addOptions(groupSelect, Object.keys(tfa.groups), "All samples");

function matchedPairs() {
  const indices = [];
  tfa.pairs.forEach((pair, index) => {
    if ((!taxonSelect.value || pair.taxon === taxonSelect.value) &&
        (!functionSelect.value || pair.function === functionSelect.value)) {
      indices.push(index);
    }
  });
  return indices;
}

function sampleValues(pairIndices) {
  const values = new Float64Array(tfa.samples.length);
  for (const pairIndex of pairIndices) {
    for (const [, sampleIndex, load] of loadsByPair[pairIndex]) {
      values[sampleIndex] += load;
    }
  }
  const grouping = tfa.groups[groupSelect.value];
  return tfa.samples.map((sample, index) => ({
    sample,
    group: grouping ? grouping[index] : "All samples",
    load: values[index]
  }));
}

function plotWidth(id) {
  const gutter = {heatmap: 170, boxplot: 90, samples: 180}[id];
  return Math.max(240, document.getElementById(id).parentElement.clientWidth - gutter);
}

function heatmapSpec(indices) {
  const top = indices
    .map(index => tfa.pairs[index])
    .sort((a, b) => b.total - a.total)
    .slice(0, Number(limitSelect.value));
  const taxaOrder = [...new Set(top.map(pair => pair.taxon))];
  const functionOrder = [...new Set(top.map(pair => pair.function))];
  return {
    $schema: schema,
    data: {values: top},
    width: Math.min(plotWidth("heatmap"), Math.max(300, functionOrder.length * 52)),
    height: Math.max(100, Math.min(560, taxaOrder.length * 29)),
    mark: {type: "rect", cursor: "pointer", stroke: "white", strokeWidth: 1},
    encoding: {
      x: {field: "function", type: "nominal", sort: functionOrder, title: "Function"},
      y: {field: "taxon", type: "nominal", sort: taxaOrder, title: "Taxon"},
      color: {field: "total", type: "quantitative", title: "Total load", scale: {scheme: "blues"}},
      tooltip: [
        {field: "taxon", type: "nominal", title: "Taxon"},
        {field: "function", type: "nominal", title: "Function"},
        {field: "total", type: "quantitative", title: "Total load", format: ".4g"}
      ]
    },
    config: {view: {stroke: null}}
  };
}

function boxplotSpec(samples) {
  return {
    $schema: schema,
    data: {values: samples},
    width: plotWidth("boxplot"),
    height: 300,
    mark: {type: "boxplot", extent: 1.5},
    encoding: {
      x: {field: "group", type: "nominal", title: groupSelect.value || "Samples"},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: {field: "group", type: "nominal", legend: null}
    },
    config: {view: {stroke: null}}
  };
}

function sampleSpec(samples) {
  const top = [...samples].sort((a, b) => b.load - a.load).slice(0, 100);
  return {
    $schema: schema,
    data: {values: top},
    width: plotWidth("samples"),
    height: 300,
    mark: {type: "bar"},
    encoding: {
      x: {field: "sample", type: "nominal", sort: "-y", title: "Sample", axis: {labels: top.length <= 30}},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: {field: "group", type: "nominal", title: groupSelect.value || "Group", legend: groupSelect.value ? {} : null},
      tooltip: [
        {field: "sample", type: "nominal", title: "Sample"},
        {field: "group", type: "nominal", title: groupSelect.value || "Group"},
        {field: "load", type: "quantitative", title: "TFA load", format: ".4g"}
      ]
    },
    config: {view: {stroke: null}}
  };
}

async function draw() {
  for (const view of activeViews) view.finalize();
  activeViews = [];
  const indices = matchedPairs();
  const samples = sampleValues(indices);
  const total = samples.reduce((sum, row) => sum + row.load, 0);
  status.textContent = `${indices.length} matching taxon–function pairs · ${samples.length} samples · total load ${total.toLocaleString()}`;
  document.getElementById("group-caption").textContent = groupSelect.value
    ? `Each box compares sample loads for one ${groupSelect.value} category. Zeros are included.`
    : "One box summarizes all samples. Add categorical metadata to compare groups.";

  if (!indices.length || !samples.length) {
    for (const id of ["heatmap", "boxplot", "samples"]) {
      document.getElementById(id).textContent = "No values to plot for this selection.";
    }
    return;
  }
  if (typeof vegaEmbed !== "function") {
    status.textContent += " · Vega could not load; an internet connection is needed to display the plots.";
    return;
  }
  try {
    const options = {actions: false, renderer: "svg"};
    const [heatmap, boxplot, samplePlot] = await Promise.all([
      vegaEmbed("#heatmap", heatmapSpec(indices), options),
      vegaEmbed("#boxplot", boxplotSpec(samples), options),
      vegaEmbed("#samples", sampleSpec(samples), options)
    ]);
    activeViews = [heatmap.view, boxplot.view, samplePlot.view];
    heatmap.view.addEventListener("click", (_, item) => {
      const pair = item && item.datum;
      if (pair && pair.taxon && pair.function) {
        taxonSelect.value = pair.taxon;
        functionSelect.value = pair.function;
        update();
      }
    });
  } catch (error) {
    status.textContent = `Could not draw the plots: ${error.message}`;
  }
}

function update() {
  renderQueue = renderQueue.then(draw).catch(error => {
    status.textContent = `Could not draw the plots: ${error.message}`;
  });
}

for (const select of [taxonSelect, functionSelect, groupSelect, limitSelect]) {
  select.addEventListener("change", update);
}
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(update, 150);
});
update();
