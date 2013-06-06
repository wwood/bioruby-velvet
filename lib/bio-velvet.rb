require 'bio'
require 'bio-logger'
Bio::Log::LoggerPlus.new('bio-velvet')
module Bio::Velvet
  module Logging
    def log
      Bio::Log::LoggerPlus['bio-velvet']
    end
  end
end

require 'bio-velvet/graph.rb'
require 'bio-velvet/runner.rb'


