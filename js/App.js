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
  T_cells: 'Goldenrod'
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
      .getElementById('bar-chart-div')
      .addEventListener("alignmentjs_wheel_event", function(e) {
        $('#bar-chart-div').scrollLeft(e.detail.x_pixel);
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
    d3.json(`/data/${patient_id}_cFEL/${patient_id}.json`, (err, json) => {
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
          <NavItem onClick={()=>this.takeSnapshot()}>Snapshot</NavItem>
          {self.patient_ids.map(patient_id => {
            return (<NavItem
              eventKey={patient_id}
              key={patient_id}
              active={patient_id == self.state.patient}
              onSelect={(eventKey)=>this.setState({patient_id: patient_id})}
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
    const { full_pixel_width, full_pixel_height, label_width } = this;
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
        "bar-chart-div"
      ]
    });
  }
  handleBarWheel(e) {
    e.preventDefault();
    this.scroll_broadcaster.handleWheel(e, 'main');
  }
  initializeBar(hyphy) {
    const data = JSON.parse(hyphy).MLE.content["0"];
    const { full_pixel_width } = this,
      { site_size } = this.props,
      site_scale = d3.scale.linear()
        .domain([1, this.number_of_sites])
        .range([site_size/2, full_pixel_width-site_size/2]),
      site_axis = d3.svg.axis()
        .scale(site_scale)
        .orient("bottom")
        .tickValues(d3.range(1, this.number_of_sites, 2));

    const plot_svg = d3.select("#bar-chart");
    plot_svg.html("");
    plot_svg.attr("width", full_pixel_width)
      .attr("height", 500);

    plot_svg
      .append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0, 475)")
      .call(site_axis);

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

    const statistic_scale = d3.scale.linear()
        .domain([y_min, y_max])
        .range([400, 0]),
      statistic_axis = d3.svg.axis()
        .scale(statistic_scale)
        .orient("left"),
      bar_axis_width = this.column_sizes[1],
      bar_axis_height = this.row_sizes[0],
      bar_axis_svg = d3.select("#bar-axis")
        .attr("width", bar_axis_width)
        .attr("height", bar_axis_height);

      bar_axis_svg.append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${bar_axis_width-1}, 75)`)
        .call(statistic_axis);
      
    bar_axis_svg.append('g')
      .attr('transform', 'translate(70, 350)')
      .append("text")
        .attr('x', 0)
        .attr('y', 0)
        .text('Magical evolutionary statistic')
        .attr('transform', 'rotate(270)');
  }
  initializeStructure(props) {
    const options = {
        width: this.column_sizes[0],
        height: this.row_sizes[0],
        antialias: true,
        quality : 'medium',
        background: '#FFF'
      },
      structure_div = document.getElementById('structure'),
      viewer = pv.Viewer(structure_div, options),
      structure = pv.io.pdb(props.structure, options),
      chain = structure.select({chain: 'A'});
    const geom = viewer.cartoon('protein', chain);
    viewer.autoZoom();
    geom.colorBy(new pv.color.ColorOp(function(atom, out, index) {
      const background_color = .8;
      out[index] = background_color;
      out[index + 1] = background_color;
      out[index + 2] = background_color;
      out[index + 3] = background_color;
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
        var index = atom.residue().num();
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

      <div id='bar-axis-svg'>
        <svg id='bar-axis' />
      </div>

      <div
        id='bar-chart-div'
        style={{overflowX: "scroll"}}
        onWheel={e => this.handleBarWheel(e)}
      >
        <svg
          id='bar-chart'
          width={this.full_pixel_width}
          height={this.row_sizes[1]}
        />
      </div>

      <div>
        <svg width={this.column_sizes[0]} height={this.row_sizes[1]}>
          <g transform="translate(60, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Monocytes} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Monocytes</text>
          </g>
          <g transform="translate(260, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.Plasma} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">Plasma</text>
          </g>
          <g transform="translate(460, 10)">
            <rect x="0" y="0" width="20" height="20" fill={legend.T_cells} />
            <text x="25" y="10" textAnchor="start" alignmentBaseline="middle">T cell</text>
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

