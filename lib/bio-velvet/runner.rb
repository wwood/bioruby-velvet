require 'files'
require 'systemu'
require 'parallel'

module Bio
  module Velvet
    class Runner
      include Bio::Velvet::Logging
      include Parallel::ProcessorCount

      # Run velveth and then velvetg, with the given kmer size. Returned
      # is a Bio::Velvet::Result class, stored in a temporary directory.
      # The temporary directory is removed upon program exit.
      #
      # The velveth_options and velvetg_options are strings to pass as arguments
      # to velveth and velvetg, respectively.
      #
      # The final options argument is used to specify bio-velvet wrapper options. Currently:
      # :output_assembly_path: a directory where the assembly takes place [default: a temporary directory]
      # :velveth_path: path to the velveth binary [default: 'velveth']
      # :velvetg_path: path to the velvetg binary [default: 'velvetg']
      # :threads: number of threads to use [default: all threads on the machine]
      def velvet(kmer_length, velveth_options_string, velvetg_options_string='', options={})
        res = velveth kmer_length, velveth_options_string, options
        velvetg res, velvetg_options_string, options
      end

      # Run velveth with the given kmer and return a Bio::Velvet::Result object
      #
      # Options:
      # :velveth_path: path to the velveth binary [default: 'velveth']
      def velveth(kmer_length, velveth_arguments, options={})
        set_num_cpus(options[:threads])
        result = Result.new
        outdir = nil
        if options[:output_assembly_path]
          log.debug "Using pre-defined assembly directory: #{options[:output_assembly_path] }"
          outdir = options[:output_assembly_path]
        else
          outdir = Files.create.root
        end
        result.result_directory = outdir

        binary = options[:velveth_path]
        binary ||= 'velveth'

        # Run velveth
        cmd = "#{binary} #{result.result_directory} #{kmer_length} #{velveth_arguments}"
        log.info "Running velveth: #{cmd}" if log.info?
        status, stdout, stderr = systemu cmd
        if status.exitstatus != 0
          raise VelvetRunnerException, "Error running velveth: #{stderr}\n#{stdout}"
        end
        result.velveth_stdout = stdout
        result.velveth_stderr = stderr

        return result
      end

      # Run velvetg, with a Bio::Velvet::Result object
      # generated with velveth, and velvetg arguments as a String (no need to specify the velvet directory, just the extra
      # arguments).
      #
      # Further options (the third argument):
      # :velvetg_path: path to the velvetg binary [default: 'velvetg']
      def velvetg(velveth_result_object, velvetg_arguments, options={})
        binary = options[:velvetg_path]
        binary ||= 'velvetg'

        cmd = "#{binary} #{velveth_result_object.result_directory} #{velvetg_arguments}"

        log.info "Running velvetg: #{cmd}" if log.info?
        status, stdout, stderr = systemu cmd
        if status.exitstatus != 0
          raise VelvetRunnerException, "Error running velvetg: #{stderr}\n#{stdout}"
        end
        velveth_result_object.velvetg_stdout = stdout
        velveth_result_object.velvetg_stderr = stderr

        return velveth_result_object
      end

      # Detect the binary version currently in use and return
      # as a String
      def binary_version
        cmd = 'velveth'
        log.info "Running velveth: #{cmd}" if log.info?
        status, stdout, stderr = systemu cmd
        if status.exitstatus != 0
          raise VelvetRunnerException, "Error running velveth: #{stderr}\n#{stdout}"
        end
        splits = stdout.split("\n")
        if splits.length > 1 and matches = splits[1].match(/^Version (.+)$/)
          return matches[1]
        else
          raise "Unable to parse the version number from running `#{cmd}', the output was: #{stdout}"
        end
      end

      # According to the manual,
      #"Velvet will the use up to OMP_NUM_THREADS+1 or OMP_THREAD_LIMIT threads."
      #
      # Argument is the number of CPUs, or nil if all CPUs on the machine are
      # to be used
      def set_num_cpus(num_cpus=nil)
        num_cpus ||= processor_count #from the parallel gem
        log.debug "Setting number of CPUs to run velvet with to #{num_cpus}."
        ENV['OMP_NUM_THREADS'] = (num_cpus - 1).to_s
        ENV['OMP_THREAD_LIMIT'] = num_cpus.to_s
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

      # Return a Bio::Velvet::Graph object built from the LastGraph file.
      # The options for parsing are as per Bio::Velvet::Graph#parse_from_file
      def last_graph(options=nil)
        Bio::Velvet::Graph.parse_from_file(last_graph_path, options)
      end
    end
  end
end
