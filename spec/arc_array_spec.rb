require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

include Bio::Velvet

describe "ArcArray" do
  it 'should push' do
    node1 = Graph::Node.new
    node1.node_id = 1
    node2 = Graph::Node.new
    node2.node_id = 2
    arc = Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2

    ary = Graph::ArcArray.new
    ary.push arc
    ary.to_a.should == [arc]

    arc2 = Graph::Arc.new
    arc2.begin_node_id = 1
    arc2.end_node_id = 2

    ary.push arc2
    ary.to_a.should == [arc,arc2]
  end

  it 'should get arcs right' do
    node1 = Graph::Node.new
    node1.node_id = 1
    node2 = Graph::Node.new
    node2.node_id = 2
    arc = Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2

    ary = Graph::ArcArray.new
    ary.push arc
    ary.get_arcs_by_node_id(1,2).should == [arc]
    ary.get_arcs_by_node_id(2,1).should == [arc]
    ary.get_arcs_by_node_id(1).should == [arc]
    ary.get_arcs_by_node_id(2).should == [arc]
    ary.get_arcs_by_node_id(3).should == []

    arc2 = Graph::Arc.new
    arc2.begin_node_id = 1
    arc2.end_node_id = 2
    ary.push arc2
    ary.get_arcs_by_node_id(1,2).should == [arc, arc2]
    ary.get_arcs_by_node_id(2,1).should == [arc, arc2]
    ary.get_arcs_by_node_id(1).should == [arc, arc2]
    ary.get_arcs_by_node_id(2).should == [arc, arc2]
    ary.get_arcs_by_node_id(3).should == []


    ary = Graph::ArcArray.new
    arc = Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 1
    ary.push arc
    ary.get_arcs_by_node_id(1).should == [arc]
    ary.get_arcs_by_node_id(1,1).should == [arc]

    arc2 = Graph::Arc.new
    arc2.begin_node_id = 1
    arc2.end_node_id = 2
    ary.push arc2
    ary.get_arcs_by_node_id(1).should == [arc, arc2]
    ary.get_arcs_by_node_id(2).should == [arc2]
  end


  it 'should length' do
    node1 = Graph::Node.new
    node1.node_id = 1
    node2 = Graph::Node.new
    node2.node_id = 2
    arc = Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2

    ary = Graph::ArcArray.new
    ary.length.should == 0
    ary.push arc
    ary.length.should == 1

    arc2 = Graph::Arc.new
    arc2.begin_node_id = 1
    arc2.end_node_id = 2

    ary.push arc2
    ary.length.should == 2
  end

  it 'should delete' do
    node1 = Graph::Node.new
    node1.node_id = 1
    node2 = Graph::Node.new
    node2.node_id = 2
    arc = Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2

    ary = Graph::ArcArray.new
    ary.push arc
    ary.to_a.should == [arc]
    ary.delete arc
    ary.to_a.should == []

    ary.push arc

    arc2 = Graph::Arc.new
    arc2.begin_node_id = 1
    arc2.end_node_id = 2
    ary.push arc2
    ary.length.should == 2
    ary.delete arc2
    ary.to_a.should == [arc]
    ary.push arc2
    ary.delete arc
    ary.to_a.should == [arc2]
  end
end
