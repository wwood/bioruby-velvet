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

Bio::Velvet::Runner.new.binary_version #=> e.g. "1.2.08"
```

By default, the ```velvet``` method passes no parameters to ```velvetg``` other than the velvet directory created by velveth. This directory is a temporary directory by default, but this can also be set. For instance, to run velvet using with a ```-cov_cutoff``` parameter in the ```velvet_dir``` directory:
```ruby
velvet_result = Bio::Velvet::Runner.new.velvet(87,
  '-short /path/to/reads.fa',
  '-cov_cutoff 3.5', 
  :output_assembly_path => 'velvet_dir')
```

The graph file can be parsed from a ```velvet_result```:
```ruby
graph = velvet_result.last_graph #=> Bio::Velvet::Graph object
```
In my experience (mostly on complex metagenomes), the graph object itself does not take as much RAM as initially expected. Most of the hard work has already been done by velvet itself, particularly if the ```-cov_cutoff``` has been set. However parsing in the graph can take many minutes or even hours if the LastGraph file is big (>500MB). The slowest part of parsing is parsing in the positions of reads i.e. using the ```-read_trkg yes``` velvet option. To speed up that process one can use e.g.
```ruby
velvet_result.last_graph(:interesting_read_ids => Set.new([1,2,3]))
``` 
To only parse read in the positions of the first 3 reads.

With a parsed graph (a ```Bio::Velvet::Graph``` object) you can interact with the graph e.g.
```ruby
graph.kmer_length #=> 87
graph.nodes #=> Bio::Velvet::Graph::NodeArray object
graph.nodes[3] #=> Bio::Velvet::Graph::Node object with node ID 3
graph.get_arcs_by_node_id(1, 3) #=> an array of arcs between nodes 1 and 3 (Bio::Velvet::Graph::Arc objects)
graph.nodes[5].noded_reads #=> array of Bio::Velvet::Graph::NodedRead objects, for read tracking
```
There is much more that can be done to interact with the graph object and its components - see the [rubydoc](http://rubydoc.info/gems/bio-velvet/Bio/Velvet/Graph).

### Parsers for `Sequences` and `CnyUnifiedSeq.names` files
With default parameters velvet generates a `Seqeunces` file, that includes read ID information and the sequences themselves.
```ruby
seqs = Bio::Velvet::Sequences.parse_from_file(File.join velvet_result.result_directory, 'Sequences')
seqs[1] => 'AAAATTGTCAGACTAGCTATCAGCATATCAGCGCGCATCTCAGACGAGCACTATC'
```
If the `-create_binary` flag is set when running `velveth`, a names file is generated that encodes the read names and IDs.
```ruby
entries = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries(
  File.join(velvet_result.result_directory, 'CnyUnifiedSeq.names'),
  ['read1','read2']
  ) #=> Hash of read name to Array of CnyUnifiedSeqNamesFileEntry objects
entries['read1'] #=> Array of CnyUnifiedSeqNamesFileEntry objects
entries['read1'][0].read_id #=> 1 (i.e. '1'.to_i)
```
When speed is required, grep can come to the rescue (at the cost of some portability)
```ruby
entries = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries_using_grep_hack(
  File.join(velvet_result.result_directory, 'CnyUnifiedSeq.names'),
  ['read1','read2']
  ) #=> same returned object as above
```
The sequences themselves are stored in a separate file when `-create_binary` is used - an interface for this is included in the [bio-velvet_underground](https://github.com/wwood/bioruby-velvet_underground) biogem.

## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/wwood/bioruby-velvet

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

This code is currently unpublished.

## Biogems.info

This Biogem is listed at [biogems.info](http://biogems.info/index.html#bio-velvet)

## Copyright

Copyright (c) 2013 Ben J Woodcroft. See LICENSE.txt for further details.

