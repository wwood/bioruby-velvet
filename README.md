# bio-velvet

[![Build Status](https://secure.travis-ci.org/wwood/bioruby-velvet.png)](http://travis-ci.org/wwood/bioruby-velvet)

```bio-velvet``` is a [biogem](biogems.info) for interacting with the [velvet](http://www.ebi.ac.uk/~zerbino/velvet/) sequence assembler. It includes both a wrapper for the velvet executable, as well as a a parser for the 'LastGraph' format files that velvet creates. This gives access to the underlying assembly graph created by velvet.

## Installation
To install ```bio-velvet``` and its rubygem dependencies:

```sh
gem install bio-velvet
```

## Usage

To run velvet with a kmer length of 87 on a set of single ended reads in ```/path/to/reads.fa```:
```ruby
require 'bio-velvet'

velvet_result = Bio::Velvet::Runner.new.velvet(87, '-short /path/to/reads.fa') #=> Bio::Velvet::Result object

contigs_file = velvet_result.contigs_path #=> path to contigs file as a String
lastgraph_file = velvet_result.last_graph_path #=> path to last graph file as a String
```

The graph file can be then parsed. In my experience (mostly on complex metagenomes), the graph object itself does not take fantastical amounts of RAM to be stored. Most of the hard work has already been done by velvet itself. However parsing in the graph can take many minutes if the LastGraph file is big (>500MB).

Given the velvet_result from above
```ruby
graph = velvet_result.last_graph #=> Bio::Velvet::Graph object
```
With this graph you can access ingteract with the graph e.g.
```ruby
graph.kmer_length #=> 87
graph.nodes #=> Bio::Velvet::Graph::NodeArray object
graph.nodes[3] #=> Bio::Velvet::Graph::Node object with node ID 3
graph.get_arcs_by_node_id(1, 3) #=> an array of arcs between nodes 1 and 3
graph.nodes[5].noded_reads #=> array of Bio::Velvet::Graph::NodedRead objects, for read tracking
```
There is obviously more that can be done to interact with the graph object - see the [rubydoc])(

## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/wwood/bioruby-velvet

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

This code is currently unpublished.

## Biogems.info

This Biogem is published at (http://biogems.info/index.html#bio-velvet)

## Copyright

Copyright (c) 2013 Ben J Woodcroft. See LICENSE.txt for further details.

