import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import './App.css';
import {
  BaseAlignment, BaseTree, fastaParser, ScrollBroadcaster, SequenceAxis,
  sortFASTAAndNewick, computeLabelWidth
} from 'alignment.js';
import { Navbar, Nav, NavDropdown, NavItem, Grid, Row, Col } from 'react-bootstrap';
import html2canvas from 'html2canvas';
import 'bootstrap/dist/css/bootstrap.css'; 
import 'bootstrap/dist/css/bootstrap-theme.css';
import 'alignment.js/lib/alignment.css';

const pv = require('bio-pv');
const d3 = require("d3");
const $ = require("jquery");
require("phylotree");


const legend = {
  Monocytes: 'red',
  Plasma: 'blue',
  T_cells: 'Goldenrod',
  Background: 'DarkGrey',
  Synonymous: 'black',
  "Monocytes vs Plasma": "purple",
  "Monocytes vs T_cells": "orange",
  "Plasma vs T_cells": "green"
};

const rgb_legend = {
  "Monocytes vs Plasma": {R: 128/255, G: 0/255, B: 128/255},
  "Monocytes vs T_cells": {R: 255/255, G: 150/255, B: 0/255},
  "Plasma vs T_cells": {R: 0/255, G: 128/255, B: 0/255}
};

const threeToOne = {
  'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
  'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
  'NAG': '*', 'FUL': '*'
};


class App extends Component {
  constructor(props) {
    super(props);
    this.state = {
      patient_id: null,
      fasta: null,
      newick: null,
      hyphy: null,
      structure: null
    };
    this.patient_ids = ['P01', 'P02', 'P13'];
  }
  componentDidMount() {
    this.fetchPatientData('P01');
    document
      .getElementById('hyphy-chart-div')
      .addEventListener("alignmentjs_wheel_event", function(e) {
        $('#hyphy-chart-div').scrollLeft(e.detail.x_pixel);
      });

  }
  shouldComponentUpdate(nextProps, nextState) {
    const different_ids = this.state.patient_id != nextState.patient_id;
    return different_ids;
  }
  componentDidUpdate(prevProps, prevState) {
    const { patient_id } = this.state;
    if(patient_id != prevState.patent) {
      this.fetchPatientData(patient_id);
    }
  }
  fetchPatientData(patient_id) {
    d3.json(`/data/${patient_id}/dashboard.json`, (err, json) => {
      this.setState({
        patient_id: patient_id,
        fasta: json.fasta,
        newick: json.newick,
        hyphy: json.hyphy,
        structure: json.structure
      });
    });
  }
  takeSnapshot() {
    html2canvas(document.getElementById('full-viz'))
      .then(function(canvas) {
        var a = document.createElement('a');
        a.href = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
        a.download = 'somefilename.png';
        a.click();
        a.remove();
      });
  }
  render() {
    const self = this;
    return (<div>
      <Navbar>
        <Navbar.Header>
          <Navbar.Brand>
            <a href="#home">ACME HyPhy Structural Viewer</a>
          </Navbar.Brand>
        </Navbar.Header>
        <Nav pullRight>
          { /*<NavItem onClick={()=>this.takeSnapshot()}>Snapshot</NavItem>*/ }

          {self.patient_ids.map(patient_id => {
            return (<NavItem
              eventKey={patient_id}
              key={patient_id}
              active={patient_id == self.state.patient_id}
              onSelect={(eventKey)=>this.fetchPatientData(patient_id)}
            >
              {patient_id}
            </NavItem>);
          })}
        </Nav>
      </Navbar>
      <Grid>
        <Row>
          <Col>
            <StructuralViz
              fasta={this.state.fasta}
              newick={this.state.newick}
              hyphy={this.state.hyphy}
              structure={this.state.structure}
              patient_id={this.state.patient}
            />
          </Col>
        </Row>
      </Grid>
    </div>);
  }
}

class StructuralViz extends Component {
  constructor(props) {
    super(props);
    this.column_sizes = [500, 200, 700];
    this.row_sizes = [500, 40, 500];
    this.initialize(props);
  }
  componentWillUpdate(nextProps) {
    this.initialize(nextProps);
  }
  componentDidUpdate() {
    this.initializeBar(this.props.hyphy);
  }
  initialize(props) {
    if(props.fasta) {
      const all_sequence_data = fastaParser(props.fasta);
      all_sequence_data.sort((a, b) => {
        if(a.header.indexOf('HXB2') > -1) return -1;
        if(b.header.indexOf('HXB2') > -1) return 1;
        if(a.header.indexOf('3JWO') > -1) return -1;
        if(b.header.indexOf('3JWO') > -1) return 1;
        return 1;      
      });

      this.reference_sequence_data = all_sequence_data.slice(0,2);
      this.patient_sequence_data = all_sequence_data.slice(2);
      this.full_pixel_height = this.props.site_size * this.patient_sequence_data.length;
      this.number_of_sites = all_sequence_data[0].seq.length;
      this.full_pixel_width = this.props.site_size * this.number_of_sites;
      this.label_width = computeLabelWidth(
        all_sequence_data,
        this.props.label_padding
      ) + 7;

      this.column_sizes[0] += this.column_sizes[1] - this.label_width;
      this.column_sizes[1] = this.label_width;

      const phylotree_size = [this.full_pixel_height, this.column_sizes[0]];
      const { phylotree } = sortFASTAAndNewick(
        this.patient_sequence_data,
        props.newick,
        phylotree_size
      );
      phylotree.style_edges(function(element, data) {
        if(data.target.Monocytes) element.style('stroke', legend.Monocytes);
        if(data.target.Plasma) element.style('stroke', legend.Plasma);
        if(data.target.T_cells) element.style('stroke', legend.T_cells);
      });
      this.phylotree = phylotree;

      this.initializeStructure(props);
      this.setScrollingEvents(props);
    }
  }
  setScrollingEvents() {
    const width = this.column_sizes[2],
      height = this.row_sizes[2];
    const { full_pixel_width, full_pixel_height } = this;
    this.scroll_broadcaster = new ScrollBroadcaster({
      width: full_pixel_width,
      height: full_pixel_height,
      x_pad: width,
      y_pad: height,
      x_pixel: this.x_pixel || 0,
      y_pixel: this.y_pixel || 0,
      bidirectional: [
        "alignmentjs-tree-div",
        "alignmentjs-labels-div",
        "reference-alignment",
        "alignmentjs-alignment",
        "hyphy-chart-div"
      ]
    });
  }
  handleHyPhyWheel(e) {
    e.preventDefault();
    this.scroll_broadcaster.handleWheel(e, 'main');
  }
  initializeBar(hyphy) {
    const { full_pixel_width } = this,
      { site_size } = this.props,
      site_scale = d3.scale.linear()
        .domain([1, this.number_of_sites])
        .range([site_size/2, full_pixel_width-site_size/2]),
      site_axis = d3.svg.axis()
        .scale(site_scale)
        .orient("bottom")
        .tickValues(d3.range(1, this.number_of_sites, 2));

    const plot_svg = d3.select("#hyphy-chart");
    plot_svg.html("");
    plot_svg.attr("width", full_pixel_width)
      .attr("height", 500);

    plot_svg
      .append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0, 475)")
      .call(site_axis);

    const y_min = 0;
    const y_max = d3.max([
      d3.max(hyphy.alpha),
      d3.max(hyphy['beta (T_cells)']),
      d3.max(hyphy['beta (Monocytes)']),
      d3.max(hyphy['beta (background)']),
      d3.max(hyphy['beta (Plasma)'])
    ]);

    const hyphy_axis_width = this.column_sizes[1],
      hyphy_axis_height = this.row_sizes[0],
      line_plot_height = 300,
      annotation_plot_height = hyphy_axis_height - line_plot_height,
      x_axis_height = 25,
      line_plot_offset = annotation_plot_height-x_axis_height,
      statistic_scale = d3.scale.linear()
        .domain([y_min, y_max])
        .range([line_plot_height, 0]),
      statistic_axis = d3.svg.axis()
        .scale(statistic_scale)
        .orient("left"),

      hyphy_axis_svg = d3.select("#hyphy-axis")
        .attr("width", hyphy_axis_width)
        .attr("height", hyphy_axis_height);
      hyphy_axis_svg.html('');

      hyphy_axis_svg.append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${hyphy_axis_width-1}, ${line_plot_offset})`)
        .call(statistic_axis);
      
    hyphy_axis_svg.append('g')
      .attr('transform', 'translate(60, 350)')
      .append("text")
        .attr('x', 0)
        .attr('y', 0)
        .text('Evolutionary rate')
        .attr('transform', 'rotate(270)');

    const categories = [
      { id: 'alpha', color: 'black' },
      { id: 'beta (background)', color: 'DarkGrey' },
      { id: 'beta (T_cells)', color: legend['T_cells'] },
      { id: 'beta (Monocytes)', color: legend['Monocytes'] },
      { id: 'beta (Plasma)', color: legend['Plasma'] }
    ];
    const plot_g = plot_svg.append('g')
      .attr("transform", `translate(20, ${line_plot_offset})`);

    categories.forEach(function(category) {
      var line = d3.svg.line()
        .x(function(d,i) { return site_scale(i); })
        .y(function(d,i) { return statistic_scale(d); });

      plot_g.append("path")
        .attr("d", line(hyphy[category.id]))
        .style("stroke", category.color)
        .style("stroke-width", "2px")
        .style("fill", "none");
    });

    const annotation_tiers = 12;
    hyphy.hxb2.forEach(function(annotated_site, i) {
        const delta = (i%(annotation_tiers-1))*annotation_plot_height/annotation_tiers,
        y1 = line_plot_offset - delta,
        height = line_plot_height + delta,
        x = site_scale(annotated_site.index);
      plot_svg.append("line")
        .attr('x1', x)
        .attr('x2', x)
        .attr('y1', y1)
        .attr('y2', y1+height)
        .style('stroke', 'LightGrey');
      plot_svg.append('circle')
        .attr('cx', x)
        .attr('cy', y1)
        .attr('r', 3)
        .attr('fill', 'black');
    });
    hyphy.hxb2.forEach(function(annotated_site, i) {
      const delta = (i%(annotation_tiers-1))*annotation_plot_height/annotation_tiers,
        y1 = line_plot_offset - delta,
        x = site_scale(annotated_site.index);
      plot_svg.append('text')
        .attr('x', x+5)
        .attr('y', y1)
        .attr('alignment-baseline', 'middle')
        .text(annotated_site.annotation);
    });
    Object.keys(rgb_legend).forEach(function(pair) {
      hyphy['P-value for ' + pair].forEach(function(site) {
        if(site.added) {
          plot_svg.append('rect')
            .attr('x', site_scale(+site.added)-site_size/2)
            .attr('y', 0)
            .attr('width', site_size)
            .attr('height', 500)
            .attr('fill', legend[pair])
            .attr('opacity', .5);
        }
      });
    });
  }
  initializeStructure(props) {
    const options = {
        width: this.column_sizes[0],
        height: this.row_sizes[0],
        antialias: true,
        quality : 'medium',
        background: '#FFF'
      },
      structure_div = document.getElementById('structure');
    structure_div.innerHTML = '';
    const viewer = pv.Viewer(structure_div, options),
      structure = pv.io.pdb(props.structure, options),
      chain = structure.select({chain: 'A'});
    const geom = viewer.cartoon('protein', chain);
    viewer.autoZoom();
    const background_color = .9;
    geom.colorBy(new pv.color.ColorOp(function(atom, out, index) {
      var r_color = background_color,
        g_color = background_color,
        b_color = background_color;
      
      Object.keys(rgb_legend).forEach(function(pair) {
        const p_values = props.hyphy['P-value for ' + pair];
        if(p_values) {
          const resnum = atom.residue().num(),
            pdb_index = resnum - 33 - (resnum > 124 ? 198 - 124 - 1: 0),
            should_highlight = p_values.map(site=>site.pdb)
            .indexOf(pdb_index) > -1;
          if(should_highlight) {
            r_color = rgb_legend[pair].R;
            g_color = rgb_legend[pair].G;
            b_color = rgb_legend[pair].B;
          }
        }
      });
      out[index] = r_color;
      out[index + 1] = g_color;
      out[index + 2] = b_color;
      out[index + 3] = 1;
    }));

    function setColorForAtom(go, atom, color){
      var view = go.structure().createEmptyView();
      view.addAtom(atom);
      go.colorBy(pv.color.uniform(color), view);
    }

    var prevPicked = null;
    structure_div.addEventListener('mousemove', function(event){
      var rect = viewer.boundingClientRect();
      var picked = viewer.pick({ x : event.clientX - rect.left,
            y : event.clientY - rect.top });
      if (prevPicked !== null && picked !== null &&
        picked.target() === prevPicked.atom){
        return;
      }
      if (prevPicked !== null){
        setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
      }
      if (picked !== null){
        var atom = picked.target();
        var residue = atom.residue();
        var index = residue.num();
        console.log(index, threeToOne[residue._name]);
        //document.getElementById('changes').innerHTML = (index+1) + changes[index]
        var color = [0,0,0,0];
        picked.node().getColorForAtom(atom, color);
        prevPicked = { atom : atom, color : color, node : picked.node() };
        setColorForAtom(picked.node(), atom, 'green');
      }
      else{
        //document.getElementById('changes').innerHTML = '&nbsp;';
        prevPicked = null;
      }
      viewer.requestRedraw();
    });
  }
  render() {
    const template_css = {
      display: "grid",
      gridTemplateColumns: this.column_sizes.join("px ") + "px",
      gridTemplateRows: this.row_sizes.join("px ") + "px"
    };
    return (<div id='full-viz' style={template_css}>
      <div id='structure'/>

      <div id='hyphy-axis-svg'>
        <svg id='hyphy-axis' />
      </div>

      <div
        id='hyphy-chart-div'
        style={{overflowX: "scroll"}}
        onWheel={e => this.handleHyPhyWheel(e)}
      >
        <svg
          id='hyphy-chart'
          width={this.full_pixel_width}
          height={this.row_sizes[1]}
        />
      </div>

      <div>
        <svg width={this.column_sizes[0]} height={this.row_sizes[1]}>
          <g transform="translate(30, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Monocytes} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Monocytes</text>
          </g>
          <g transform="translate(140, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Plasma} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Plasma</text>
          </g>
          <g transform="translate(250, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.T_cells} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">T cell</text>
          </g>
          <g transform="translate(360, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Background} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Background</text>
          </g>
          <g transform="translate(470, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Synonymous} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Synonymous</text>
          </g>
        </svg>
      </div>

      <SequenceAxis
        width={this.column_sizes[1]}
        height={this.row_sizes[1]}
        sequence_data={this.reference_sequence_data}
        id="reference"
      />

      <BaseAlignment
        width={this.column_sizes[2]}
        height={this.row_sizes[1]}
        sequence_data={this.reference_sequence_data}
        id='reference'
        amino_acid
        disableVerticalScrolling
        scroll_broadcaster={this.scroll_broadcaster}
      />

      <BaseTree
        phylotree={this.phylotree}
        scroll_broadcaster={this.scroll_broadcaster}
      />

      <SequenceAxis
        width={this.column_sizes[1]}
        height={this.row_sizes[2]}
        sequence_data={this.patient_sequence_data}
        scroll_broadcaster={this.scroll_broadcaster}
      />

      <BaseAlignment
        width={this.column_sizes[2]}
        height={this.row_sizes[2]}
        amino_acid
        sequence_data={this.patient_sequence_data}
        scroll_broadcaster={this.scroll_broadcaster}
      />
    </div>);
  }
}

StructuralViz.defaultProps = {
  site_size: 20,
  label_padding: 5
}

class OldStructuralViz extends Component {
  componentDidUpdate() {
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
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)

