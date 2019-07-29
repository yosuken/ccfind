#
#  ccfind.rake - tool for finding circular complete sequences by detecting terminal redundancy of each sequence
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#

# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, queue, nthreads, mem, wtime, ncpus|
	# [TODO] queue validation
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if queue != ""
			raise("`--queue #{queue}': invalid queue") unless %w|JP1 JP4 JP10 cdb|.include?(queue)
			sh "qsubarraywww -q #{queue} -l ncpus=#{nthreads} -l mem=#{mem}gb -l walltime=#{wtime} #{fin}"
		elsif ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/
			sh "parallel --jobs #{ncpus} <#{fin}"
		else
			sh "sh #{fin}"
		end
	}
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

CheckVersion = lambda do |commands|
	commands.each{ |command|
		str = case command
					when "makeblastdb"
						%|makeblastdb -version 2>&1|
					when "blastn"
						%|blastn -version 2>&1|
					when "ssearch"
						%{ssearch 2>&1 |head -n 7 |tail -n 3}
					when "ruby"
						%|ruby --version 2>&1|
					when "parallel"
						%{LANG=C parallel --version 2>&1 |head -n 1}
					when "prodigal"
						%{prodigal 2>&1 |head -n 2 |tail -n 1}
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "$ #{str}"
		### run
		puts `#{str}`
		### flush
		$stdout.flush
	}
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	### define tasks
  tasks  = %w|01-1.trim_fasta 01-2.makeblastdb 02.blastn 03.parse_blastn|
  tasks += %w|04-1.prepare_ssearch 04-2.ssearch 04-3.parse_ssearch 04-4.aln_extend 04-5.make_circ_list_and_fasta|
  tasks += %w|04-6.prodigal_for_circ 04-7.move_start_position_for_circ| ## optional

  ### constants
  Fin    = ENV["fin"]
  Size   = ENV["size"].to_i
  Len    = ENV["len"].to_i
  Idt    = ENV["idt"].to_i

  dir    = ENV["dir"]
  Sdir   = "#{dir}/tmp/S"
  Edir   = "#{dir}/tmp/E"
  Odir   = "#{dir}/tmp/O"
  Rdir   = "#{dir}/result"
  Ridir  = "#{dir}/result/intermediate"
  Jdir   = "#{dir}/batch"

  Sfa    = "#{Sdir}/S#{Size}.fasta"
  Efa    = "#{Edir}/E#{Size}.fasta"
  Otab   = "#{Odir}/blastn.out"
  Olog   = "#{Odir}/blastn.log"

  Dbsize = 10_000_000 ## blastn -dbsize
  MaxTS  = 10_000_000 ## baastn -max_target_seqs
  Evalue = ENV["evalue"]||"10"

	Qname  = ENV["queue"]||""
	Wtime  = ENV["wtime"]||"24:00:00"
	Ncpus  = ENV["ncpus"]||""    
	Nthreads = "1".to_i              
	Mem      = Nthreads * 12

	### check version
	commands  = %w|makeblastdb blastn ssearch ruby|
	commands += %w|parallel| if Ncpus != ""
	CheckVersion.call(commands)

	### run
	NumStep  = tasks.size
	tasks.each.with_index(1){ |task, idx|
		Rake::Task[task].invoke(idx)
	}
end
# }}} default


# {{{ tasks
desc "01-1.trim_fasta"
task "01-1.trim_fasta", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	# [!!!] sequences are skipped if length is < (2 x Size).
  [Odir, Sdir, Edir, Rdir, Ridir].each{ |_dir| mkdir_p _dir unless File.directory?(_dir) }

  outS = []
  outE = []
	open("#{Ridir}/01.skipped.list", "w"){ |fskip|
		IO.read(Fin).split(/^>/)[1..-1].each{ |ent|
			lab, *seq = ent.split(/\n/)
			seq = seq*""
			if seq.size < 2 * Size # length
				fskip.puts [lab, "#{seq.size} nt (< #{2 * Size} nt)"]*"\t"
				next 
			end
      outS << [">"+lab, seq[0, Size]]*"\n"
      outE << [">"+lab, seq[-Size..-1]]*"\n"
		}
	}
  open(Sfa, "w"){ |foutS| foutS.puts outS }
  open(Efa, "w"){ |foutE| foutE.puts outE }
end
desc "01-2.makeblastdb"
task "01-2.makeblastdb", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  command = "makeblastdb -dbtype nucl -in #{Efa} -out #{Efa} -title #{File.basename(Efa)} 2>#{Efa}.makeblastdb.log"
  sh command
end
desc "02.blastn"
task "02.blastn", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  command = "blastn -dust no -word_size 12 -strand plus -dbsize #{Dbsize} -max_target_seqs #{MaxTS} -num_threads #{Ncpus == "" ? 1 : Ncpus} -evalue #{Evalue} -outfmt 6 -db #{Efa} -query #{Sfa} -out #{Otab} 2>#{Olog}"
  sh command
end
desc "03.parse_blastn"
task "03.parse_blastn", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  open("#{Ridir}/03.blastn.hit.out", "w"){ |fout|
    IO.readlines(Otab).each{ |l| 
      next if l =~ /^#/
      a = l.chomp.split("\t")
      next if a[0] != a[1]
      fout.puts l 
    }
  }
end
desc "04-1.prepare_ssearch"
task "04-1.prepare_ssearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  jdir = "#{Jdir}/04-1"; mkdir_p jdir

  sdir = "#{Sdir}/each"
  edir = "#{Edir}/each"
  odir = "#{Odir}/each"
	[odir, sdir, edir].each{ |_dir| mkdir_p _dir unless File.directory?(_dir) }

  labs = {} ### store BLASTn-positive sequences
  IO.readlines("#{Ridir}/03.blastn.hit.out").each{ |l|
    lab = l.chomp.split("\t")[0]
    labs[lab] = 1
  }

  [Sfa, Efa].zip([sdir, edir]){ |fa, dir|
    IO.read(fa).split(/^>/)[1..-1].each{ |ent|
      lab, *seq = ent.split("\n")
      next unless labs[lab] ### skip if no hit by BLASTn
      open("#{dir}/#{lab}.fasta", "w"){ |fw|
        fw.puts [">"+lab, seq]
      }
    }
  }

  outs = []
  labs.each_key{ |lab|
    sfile = "#{sdir}/#{lab}.fasta"
    efile = "#{edir}/#{lab}.fasta"
    ofile = "#{odir}/#{lab}.ssearch.out"
		command = "ssearch36 -3 -n -m 8 -T 1 #{sfile} #{efile} >#{ofile}"
    outs << command
  }

	WriteBatch.call(outs, jdir, t)
end
desc "04-2.ssearch"
task "04-2.ssearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  jdir = "#{Jdir}/04-1" ## output of 04-1
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
desc "04-3.parse_ssearch"
task "04-3.parse_ssearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  fout1 = open("#{Ridir}/03.ssearch.all.out", "w")
  fout2 = open("#{Ridir}/03.ssearch.aln.out", "w")
  fout3 = open("#{Ridir}/03.ssearch.aln.idt80.out", "w")

	Dir["#{Odir}/each/*.ssearch.out"].each{ |fin|
		IO.readlines(fin).each{ |l| 
			next if l =~ /^#/
			fout1.puts l 
			a = l.chomp.split("\t")
			next if a[3].to_i < Len # aln_len >= Len
			fout2.puts l 
			next if a[2].to_f < 80  # %idt >= 80
			fout3.puts l 
		}
	}
	[fout1, fout2, fout3].each{ |_fout| _fout.close }
end
desc "04-4.aln_extend"
task "04-4.aln_extend", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  fin  = "#{Ridir}/03.ssearch.aln.idt80.out"
  fout = "#{Ridir}/04.ssearch.aln.idt80.extend.out"
  open(fout, "w"){ |fw|
    IO.readlines(fin).each{ |l|
      a = l.chomp.split("\t")
      aln_len, mismatch, gap, q_start, q_stop, s_start, s_stop = a.values_at(3..9).map(&:to_i)

      # calc offset
      start_offset = q_start - 1
      stop_offset  = Size - s_stop

      # fix alignment
      s_start  -= start_offset
      q_start  -= start_offset
      s_stop   += stop_offset
      q_stop   += stop_offset
      aln_len  += start_offset + stop_offset
      mismatch += start_offset + stop_offset

      # validation1
      ($stderr.puts [1, a]*"\t"; next) if s_start < 1
      # ($stderr.puts [2, a]*"\t"; next) if q_stop  > Size

      idt = (aln_len - mismatch - gap).to_f / aln_len * 100
      idt = "%.2f" % idt
      fw.puts [a[0..1], idt, aln_len, mismatch, gap, q_start, q_stop, s_start, s_stop, a[10..-1]]*"\t"
    }
  }
end
desc "04-5.make_circ_list_and_fasta"
task "04-5.make_circ_list_and_fasta", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  fin0  = "#{Ridir}/04.ssearch.aln.idt80.extend.out"
	fout0 = open("#{Rdir}/circ.detected.list", "w")
	fout1 = open("#{Rdir}/circ.fasta", "w")      # fasta of circular sequence
  fout2 = open("#{Rdir}/circ.noTR.fasta", "w") # fasta of circular sequence (Terminal Redundancy at the end of sequence is trimmed)

	lab2remove_len = {}
	IO.readlines(fin0).select{ |l| l.chomp.split(/\t/)[2].to_f >= Idt }.each{ |l|
		fout0.puts l
		a = l.chomp.split(/\t/)
		lab = a[0]
		s_start, s_stop = a.values_at(8..9).map(&:to_i)
		lab2remove_len[lab] = s_stop - s_start + 1
	}

	IO.read(Fin).split(/^>/)[1..-1].each{ |ent|
		lab, *seq = ent.split("\n")

    ### skip non-circular
    next unless remove_len = lab2remove_len[lab]

    seq  = seq.join("").gsub(/\s+/, "")
    _seq = seq[0, seq.size - remove_len] # remove the last 'remove_len' seq

		fout1.puts [">"+lab,  seq]
		fout2.puts [">"+lab, _seq]
	}

  [fout0, fout1, fout2].each{ |fout| fout.close }
end
desc "04-6.prodigal_for_circ"
task "04-6.prodigal_for_circ", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  fin0  = "#{Rdir}/circ.fasta"
  fout  = "#{Ridir}/05.circ.fasta.gff"

  sh "prodigal -p meta -i #{fin0} -f gff -o #{fout}"
end
desc "04-7.move_start_position_for_circ"
task "04-7.move_start_position_for_circ", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  ### move start position to the rightmost nucleotide in the rightmost intergenic region
  fin0  = "#{Rdir}/circ.noTR.fasta"
  fin1  = "#{Rdir}/circ.fasta"
  fin   = "#{Ridir}/05.circ.fasta.gff"
  fout  = "#{Rdir}/circ.noTR.cPerm.fasta"

	script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  sh "ruby #{script} #{fin0} #{fin1} #{fin} #{fout}"
end
# }}} tasks
