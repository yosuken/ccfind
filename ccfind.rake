#
#  ccfind.rake - circular complete sequences finder
#
#    Copyright: 2017- (C) Yosuke Nishimura (ynishimura@aori.u-tokyo.ac.jp)
#    License: MIT license
#

# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
  mkdir_p(jdir, { :verbose => nil }) unless File.directory?(jdir)
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, ncpus|
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/

      ldir = "#{File.dirname(fin)}/log_parallel"
      mkdir_p ldir unless File.directory?(ldir)
      sh "parallel --jobs #{ncpus} --joblog #{ldir}/#{File.basename(fin)} <#{fin}"
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
					when "ssearch36"
						%{ssearch36 2>&1 |head -n 7 |tail -n 3}
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
  ### constants
  Fin      = ENV["fin"]
  Size     = ENV["size"].to_i
  Len      = ENV["len"].to_i
  Idt      = ENV["idt"].to_i

	Ncpus    = ENV["ncpus"]||""    

  Tdir     = "tmp"       ## tmpdir
  Idir     = "#{Tdir}/I" ## input
  Sdir     = "#{Tdir}/S" ## start region
  Edir     = "#{Tdir}/E" ## end region
  Odir     = "#{Tdir}/O" ## output
  Pdir     = "#{Tdir}/P" ## prodigal
  Rdir     = "result"
  Ridir    = "result/intermediate"
  Jdir     = "batch"

  Fa       = "#{Idir}/in.fasta"

  Nbin     = 1000       ## number of sequences to split (faster computation for big input)
  Dbsize   = 10_000_000 ## blastn -dbsize
  MaxTS    = 10_000_000 ## baastn -max_target_seqs
  Evalue   = 10         ## for blastn screening

  Keeptmp  = ENV["keeptmp"]||""

	### define tasks
  tasks  = %w|01-0.prepare_fasta 01-1.trim_fasta 01-2.makeblastdb 02.blastn 03.parse_blastn|
  tasks += %w|04-1.prepare_ssearch 04-2.ssearch 04-3.parse_ssearch 04-4.aln_extend 04-5.make_circ_list_and_fasta|
  tasks += %w|04-6.prodigal_for_circ 04-7.move_start_position_for_circ|
  tasks += %w|04-8.remove_tmpdir| if Keeptmp != "true"

	### check version
	commands  = %w|makeblastdb blastn ssearch36 ruby|
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
desc "01-0.prepare_fasta"
task "01-0.prepare_fasta", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	### [1] sequences are skipped if length is < (2 x Size).
	### [2] clean comment line (only ids)
	### [3] id ducplication check
  [Idir, Sdir, Edir, Odir, Pdir, Rdir, Ridir].each{ |_dir| mkdir_p(_dir, { :verbose => nil }) unless File.directory?(_dir) }

  ids = {}
  skip_count = 0
  open(Fa, "w"){ |fw|
    open("#{Rdir}/too_short_seq.list", "w"){ |fskip|
      IO.read(Fin).split(/^>/)[1..-1].each{ |ent|
        lab, *seq = ent.split(/\n/)

        ### [1] sequences are skipped if length is < (2 x Size).
        seq = seq*""
        if seq.size < 2 * Size # length
          fskip.puts ["#{lab}", "#{seq.size} nt (< #{2 * Size} nt)"]*"\t"
          skip_count += 1
          next 
        end

        ### [2] clean comment line (only ids)
        id = lab.split(/\s+/)[0]

        ### [3] id ducplication check
        raise("sequence name is not unique (#{id}). Aborting...") if ids[id]
        ids[id] = 1

        ### make output
        fw.puts [">"+id, seq]
      }
    }
  }
  $stderr.puts "[35m[Warning][0m #{skip_count} sequence(s) are found too short. These sequences are excluded from analysis." if skip_count > 0
end
desc "01-1.trim_fasta"
task "01-1.trim_fasta", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  outS = []
  outE = []
  IO.read(Fa).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split(/\n/)
    seq = seq*""
    outS << [">"+lab, seq[0, Size]]*"\n"
    outE << [">"+lab, seq[-Size..-1]]*"\n"
  }

  outS.each_slice(Nbin).with_index(1){ |ents, idx|
    sdir = "#{Sdir}/split#{Size}/#{idx}"
    mkdir_p(sdir, { :verbose => nil })
    open("#{sdir}/S#{Size}.fasta", "w"){ |foutS| foutS.puts ents }
  }
  outE.each_slice(Nbin).with_index(1){ |ents, idx|
    edir = "#{Edir}/split#{Size}/#{idx}"
    mkdir_p(edir, { :verbose => nil })
    open("#{edir}/E#{Size}.fasta", "w"){ |foutE| foutE.puts ents }
  }
end
desc "01-2.makeblastdb"
task "01-2.makeblastdb", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  jdir = "#{Jdir}/01-2"

  outs = []
  Dir["#{Edir}/split#{Size}/*"].sort_by{ |i| File.basename(i).to_i }.each{ |edir|
    efa = "#{edir}/E#{Size}.fasta"
    outs << "makeblastdb -dbtype nucl -in #{efa} -out #{efa} -title #{File.basename(efa)} >#{efa}.makeblastdb.log 2>&1"
  }

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Ncpus)
end
desc "02.blastn"
task "02.blastn", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  jdir = "#{Jdir}/02"

  outs = []
  Dir["#{Edir}/split#{Size}/*"].map{ |i| File.basename(i).to_i }.sort.each{ |idx|
    efa  = "#{Edir}/split#{Size}/#{idx}/E#{Size}.fasta"
    sfa  = "#{Sdir}/split#{Size}/#{idx}/S#{Size}.fasta"
    odir = "#{Odir}/split#{Size}/#{idx}"
    mkdir_p(odir, { :verbose => nil })
    outs << "blastn -dust no -word_size 12 -strand plus -dbsize #{Dbsize} -max_target_seqs #{MaxTS} -num_threads 1 -evalue #{Evalue} -outfmt 6 -db #{efa} -query #{sfa} -out #{odir}/blastn.out 2>#{odir}/blastn.log"
  }

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Ncpus)
end
desc "03.parse_blastn"
task "03.parse_blastn", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  open("#{Ridir}/03.blastn.hit.out", "w"){ |fout|
    Dir["#{Odir}/split#{Size}/*/blastn.out"].sort_by{ |i| i.split("/")[-2].to_i }.each{ |ftab|
      IO.readlines(ftab).each{ |l| 
        next if l =~ /^#/
        a = l.chomp.split("\t")
        next if a[0] != a[1]
        fout.puts l 
      }
    }
  }
end
desc "04-1.prepare_ssearch"
task "04-1.prepare_ssearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  jdir = "#{Jdir}/04-1"
  sdir = "#{Sdir}/each"
  edir = "#{Edir}/each"
  odir = "#{Odir}/each"
	[odir, sdir, edir].each{ |_dir| mkdir_p(_dir, { :verbose => nil }) unless File.directory?(_dir) }

  labs = {} ### store BLASTn-positive sequences
  IO.readlines("#{Ridir}/03.blastn.hit.out").each{ |l|
    lab = l.chomp.split("\t")[0]
    labs[lab] = 1
  }

  sfas = Dir["#{Sdir}/split#{Size}/*/S#{Size}.fasta"].sort_by{ |i| i.split("/")[-2].to_i }
  efas = Dir["#{Edir}/split#{Size}/*/E#{Size}.fasta"].sort_by{ |i| i.split("/")[-2].to_i }
  [sfas, efas].zip([sdir, edir]){ |fas, dir|
    fas.each{ |fa|
      IO.read(fa).split(/^>/)[1..-1].each{ |ent|
        lab, *seq = ent.split("\n")
        next unless labs[lab] ### skip if no hit by BLASTn
        open("#{dir}/#{lab}.fasta", "w"){ |fw|
          fw.puts [">"+lab, seq]
        }
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
	RunBatch.call(jdir, Ncpus)
end
desc "04-3.parse_ssearch"
task "04-3.parse_ssearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  fout1 = open("#{Ridir}/03.ssearch.all.out", "w")
  fout2 = open("#{Ridir}/03.ssearch.aln.out", "w")
  fout3 = open("#{Ridir}/03.ssearch.aln.idt80.out", "w")

  IO.read(Fa).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    Dir["#{Odir}/each/#{lab}.ssearch.out"].each{ |fin|
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
      ($stderr.puts [2, a]*"\t"; next) if q_stop  > Size

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

	IO.read(Fa).split(/^>/)[1..-1].each{ |ent|
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
  n      = Ncpus != "" ? Ncpus.to_i : 1
  entset = Array.new(n, "")

  if File.zero?("#{Rdir}/circ.detected.list")
    $stderr.puts "SKIP 04-6 because circular sequence is not detected."
  else
    fin0   = "#{Rdir}/circ.fasta"
    IO.read(fin0).split(/^>/)[1..-1].each.with_index{ |ent, idx|
      entset[idx % n] += ">#{ent}"
    }

    outs = []
    entset.each.with_index(1){ |ents, idx|
      fin   = "#{Pdir}/circ.#{idx}.fasta"
      open(fin, "w"){ |fw| fw.puts ents }
      fout  = "#{Pdir}/05.circ.#{idx}.fasta.gff"
      flog  = "#{Pdir}/prodigal.#{idx}.log"
      outs << "prodigal -p meta -i #{fin} -f gff -o #{fout} >#{flog} 2>&1"
    }

    jdir = "#{Jdir}/04-6"
    WriteBatch.call(outs, jdir, t)
    RunBatch.call(jdir, Ncpus)
  end
end
desc "04-7.move_start_position_for_circ"
task "04-7.move_start_position_for_circ", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  if File.zero?("#{Rdir}/circ.detected.list")
    $stderr.puts "SKIP 04-7 because circular sequence is not detected."
  else
    ### move start position to the last nucleotide in the last intergenic region
    fin0  = "#{Rdir}/circ.noTR.fasta"
    fin1  = "#{Rdir}/circ.fasta"
    fgffs = "#{Pdir}/05.circ.*.fasta.gff"
    fgff  = "#{Ridir}/05.circ.fasta.gff"
    fout  = "#{Rdir}/circ.noTR.cPerm.fasta"

    script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

    sh %|ruby #{script} #{fin0} #{fin1} '#{fgffs}' #{fgff} #{fout}|
  end
end
desc "04-8.remove_tmpdir"
task "04-8.remove_tmpdir", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  ### remove tmpdir
  rm_rf(Tdir, { :verbose => nil })
end
# }}} tasks
