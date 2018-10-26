import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import './App.css';
import { BaseAlignment, fastaParser, ScrollBroadcaster } from 'alignment.js';
import html2canvas from 'html2canvas';
const pv = require('bio-pv');
const d3 = require("d3");
require("phylotree");


class App extends Component {
  constructor(props) {
    super(props);
    this.state = { json: null };
    this.sequence_data = null;
  }
  componentDidMount() {
    d3.json("/data/P01.json", (err, json) => {
      this.sequence_data = fastaParser(json.fasta);
      this.number_of_sequences = this.sequence_data.length;
      this.site_size = 20;
      this.tree_size = this.number_of_sequences * this.site_size;
      this.tree = d3.layout
        .phylotree()
        .options({
          "left-right-spacing": "fit-to-size",
          "top-bottom-spacing": "fit-to-size",
          "show-scale": false,
          "align-tips": true,
          "show-labels": false,
          selectable: false
        })
        .size([this.tree_size, 500])
        .node_circle_size(0);
      this.parsed = d3.layout.newick_parser(json.newick);
      this.tree(this.parsed);

      var i = 0;
      this.tree.traverse_and_compute(function(n) {
        var d = 1;
        if (!n.name) {
          n.name = "Node" + i++;
        }
        if (n.children && n.children.length) {
          d += d3.max(n.children, function(d) {
            return d["count_depth"];
          });
        }
        n["count_depth"] = d;
      });

      this.tree.resort_children(function(a, b) {
        return a["count_depth"] - b["count_depth"];
      }, null, null, true);

      const ordered_leaf_names = this.tree
        .get_nodes(true)
        .filter(d3.layout.phylotree.is_leafnode)
        .map(d => d.name.split('_')[0]);

      this.sequence_data.sort((a, b) => {
        const a_index = ordered_leaf_names.indexOf(a.header),
          b_index = ordered_leaf_names.indexOf(b.header);
        return a_index - b_index;
      });
      this.setState({json: json});
    });
  }
  componentDidUpdate() {
    const options = {
        width: 500,
        height: 500,
        antialias: true,
        quality : 'medium',
        background: '#FFF'
      },
      structure_div = document.getElementById('structure'),
      viewer = pv.Viewer(structure_div, options),
      structure = pv.io.pdb(this.state.json.structure, options),
      chain = structure.select({chain: 'A'});
    viewer.cartoon('protein', chain);
    viewer.autoZoom();

    this.tree.svg(d3.select("#tree")).layout();

    const number_of_sites = this.sequence_data[0].seq.length,
        full_alignment_width = this.site_size * number_of_sites;
    const scroll_broadcaster = new ScrollBroadcaster(
      { width: full_alignment_width, height: this.tree_size },
      { width: 700, height: 500},
      { x_pixel: 0, y_pixel: 0 },
      [
        "tree-div",
        "alignmentjs-alignment",
        "plot-div"
      ]
    );
    
    scroll_broadcaster.setListeners();

    $("#tree-div").off("wheel");
    $("#tree-div").on("wheel", function(e) {
      const e_mock = {
        originalEvent: {
          deltaX: 0,
          deltaY: e.originalEvent.deltaY
        }
      };
      scroll_broadcaster.handleWheel(e_mock, "tree");
    });

    $("#alignmentjs-alignment").off("wheel");
    $("#alignmentjs-alignment").on("wheel", function(e) {
      e.preventDefault();
      scroll_broadcaster.handleWheel(e, "alignment");
    });

    document
      .getElementById("tree-div")
      .addEventListener("alignmentjs_wheel_event", function(e) {
        $("#tree-div").scrollTop(e.detail.y_pixel);
      });

    document
      .getElementById("plot-div")
      .addEventListener("alignmentjs_wheel_event", function(e) {
        $("#plot-div").scrollLeft(e.detail.x_pixel);
      });

    const x_scale = d3.scale.linear()
      .domain([1, number_of_sites])
      .range([this.site_size/2, full_alignment_width-this.site_size/2]);

    var plot_svg = d3.select("#plot");
    plot_svg.html("");
    plot_svg.attr("width", full_alignment_width)
      .attr("height", 500);

    var alignment_axis = d3.svg.axis()
      .orient("bottom")
      .scale(x_scale)
      .tickValues(d3.range(1, number_of_sites, 2));

    plot_svg
      .append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0, 475)")
      .call(alignment_axis);

    const data = JSON.parse(this.state.json.hyphy).MLE.content["0"];
    const y_min = d3.min(data.map(row=>{
      return d3.min([
        (row[1]-row[0])/row[9],
        (row[2]-row[0])/row[9],
        (row[3]-row[0])/row[9],
        (row[4]-row[0])/row[9]
      ])
    }));
    const y_max = d3.max(data.map(row=>{
      return d3.max([
        (row[1]-row[0])/row[9],
        (row[2]-row[0])/row[9],
        (row[3]-row[0])/row[9],
        (row[4]-row[0])/row[9]
      ])
    }));
    const y_scale = d3.scale.linear()
      .domain([y_min, y_max])
      .range([0, 400]);
    const tcell_line = d3.svg.line()
      .x((d,i)=>x_scale(i+1))
      .y(d=>y_scale((d[1]-d[0])/(d[9]+.001)));
    plot_svg.append("path")
      .attr("class", "line")
      .attr("d", tcell_line(data))
      .style("fill", "none")
      .style("stroke-width", 2)
      .style("stroke", "#000");

    // For html2canvas
    d3.selectAll('path')
      .attr('fill', 'none')
      .attr('stroke', '#999')
      .attr('stroke-width', '2px');
    d3.selectAll('.axis text')
      .attr('font', 'san-serif');
    d3.selectAll('.axis path')
      .attr('fill', 'none')
      .attr('stroke', '#000')
      .attr('shape-rendering', 'crispEdges');
    d3.selectAll('.axis line')
      .attr('fill', 'none')
      .attr('stroke', '#000')
      .attr('shape-rendering', 'crispEdges');
    d3.selectAll('.axis')
      .attr('shape-rendering', 'crispEdges')
      .attr('font', '10 px sans-serif');
    d3.selectAll('.branch-tracer')
      .attr('stroke', '#bbb')
      .attr('stroke-dasharray', '3,4')
      .attr('stroke-width', '1px');

  }
  takeSnapshot() {
    html2canvas(document.getElementById('App'))
      .then(function(canvas) {
        var a = document.createElement('a');
        a.href = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
        a.download = 'somefilename.png';
        a.click();
        a.remove();
      });
  }
  render() {
    return (<div>
      <button onClick={()=>this.takeSnapshot()}>Snapshot</button>
      <div id="App">
        <div id="structure" />
        <div id="plot-div">
          <svg id="plot" />
        </div>
        <div id="tree-div">
          <svg id="tree" />
        </div>
        <BaseAlignment
          amino_acid
          width={700}
          height={500}
          sequence_data={this.sequence_data}
        />
      </div>
    </div>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)

