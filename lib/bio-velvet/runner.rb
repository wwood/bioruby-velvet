require 'files'
require 'systemu'

module Bio
  module Velvet
    class Runner
      include Bio::Velvet::Logging

      # Run velveth and then velvetg, with the given kmer size. Returned
      # is a Bio::Velvet::Result class, stored in a temporary directory.
      # The temporary directory is removed upon program exit.
      #
      # The velveth_options and velvetg_options are strings to pass as arguments
      # to velveth and velvetg, respectively.
      def velvet(kmer_length, velveth_options_string, velvetg_options_string='')
        res = velveth kmer_length, velveth_options_string
        velvetg res, velvetg_options_string
      end

      def velveth(kmer_length, velveth_arguments)
        result = Result.new
        outdir = Files.create.root
        result.result_directory = outdir

        # Run velveth
        cmd = "velveth #{result.result_directory} #{kmer_length} #{velveth_arguments}"
        log.info "Running velveth: #{cmd}" if log.info?
        status, stdout, stderr = systemu cmd
        if status.exitstatus != 0
          raise VelvetRunnerException, "Error running velveth: #{stderr}"
        end
        result.velveth_stdout = stdout
        result.velveth_stderr = stderr

        return result
      end

      # Run velvetg, with a Bio::Velvet::Result object
      # generated with velveth
      def velvetg(velveth_result_object, velvetg_arguments)
        cmd = "velvetg #{velveth_result_object.result_directory} #{velvetg_arguments}"
        log.info "Running velvetg: #{cmd}" if log.info?
        status, stdout, stderr = systemu cmd
        if status.exitstatus != 0
          raise VelvetRunnerException, "Error running velvetg: #{stderr}"
        end
        velveth_result_object.velvetg_stdout = stdout
        velveth_result_object.velvetg_stderr = stderr

        return velveth_result_object
      end
    end

    class VelvetRunnerException < Exception; end

    class Result
      attr_accessor :velveth_stdout, :velveth_stderr
      attr_accessor :velvetg_stdout, :velvetg_stderr
      attr_accessor :result_directory

      # Path to the LastGraph output from velvetg
      def last_graph_path
        File.join result_directory, 'LastGraph'
      end

      # Path to the contigs.fa output from velvetg
      def contigs_path
        File.join result_directory, 'contigs.fa'
      end

      # Path to the stats.txt output from velvetg
      def stats_path
        File.join result_directory, 'stats.txt'
      end

      # Return a Bio::Velvet::Graph object built from the LastGraph file
      def last_graph
        Bio::Velvet::Graph.parse_from_file(last_graph_path)
      end
    end
  end
end
