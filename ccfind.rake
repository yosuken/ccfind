#
#  ccfind.rake - tool for finding circular complete sequences by detecting terminal redundancy of each sequence
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#

desc "01.make_term_fasta"
task "01.make_term_fasta" do |t|
	# [!!!] sequences are skipped if length is < (2 x Size).
	dir  = ENV["dir"]
	fin  = ENV["fin"]
	Size = ENV["Size"].to_i
	sdir = "#{dir}/node/S"
	edir = "#{dir}/node/E"
	odir = "#{dir}/node/O"
	[odir, sdir, edir].each{ |_dir| mkdir_p _dir unless File.directory?(_dir) }

	ldir = "#{dir}/result"
	mkdir_p ldir unless File.directory?(ldir)

	open("#{ldir}/01.skipped.list", "w"){ |fskip|
		IO.read(fin).split(/^>/)[1..-1].each{ |ent|
			label, *seq = ent.split(/\n/)
			seq = seq*""
			if seq.size < 2 * Size # length
				fskip.puts [label, "#{seq.size} nt (< #{2 * Size} nt)"]*"\t"
				next 
			end
			open("#{sdir}/#{label}.fasta", "w"){ |fout1|
				fout1.puts [">"+label, seq[0, Size]]
			}
			open("#{edir}/#{label}.fasta", "w"){ |fout2|
				fout2.puts [">"+label, seq[-Size..-1]]
			}
		}
	}
end
desc "02.ssearch"
task "02.ssearch" do |t|
	dir  = ENV["dir"]
	sdir = "#{dir}/node/S"

	Dir["#{sdir}/*.fasta"].each{ |sfile|
		node  = File.basename(sfile).split(".")[0..-2]*"."
		efile = sfile.gsub(%r{/S/}, "/E/")
		odir  = "#{dir}/node/O"
		ofile = "#{odir}/#{node}.ssearch.StoE.out"
		command = "ssearch36 -3 -n -m 8 -T 1 #{sfile} #{efile} >#{ofile}"
		$stderr.puts command

		IO.popen(command, :err => [:child, :out]){ |pipe|
			pipe.each_line{ |l|
				$stderr.puts l
			}
		}
		raise if $?.exitstatus != 0
	}
end
desc "03.parse_ssearch"
task "03.parse_ssearch" do |t|
	dir  = ENV["dir"]
	Len  = ENV["Len"].to_i
	odir = "#{dir}/result"

	fout1 = open("#{odir}/03.all.out", "w")
	fout2 = open("#{odir}/03.aln.out", "w")
	fout3 = open("#{odir}/03.aln.idt80.out", "w")
	# fout4 = open("#{odir}/03.aln.idt90.out", "w")
	Dir["#{dir}/node/O/*.ssearch.StoE.out"].each{ |fin|
		IO.readlines(fin).each{ |l| 
			next if l =~ /^#/
			fout1.puts l 
			a = l.chomp.split("\t")
			next if a[3].to_i < Len # aln_len >= Len
			fout2.puts l 
			next if a[2].to_f < 80  # %idt >= 80
			fout3.puts l 
			# next if a[2].to_f < 90  # %idt >= 90
			# fout4.puts l 
		}
	}
	# [fout1, fout2, fout3, fout4].each{ |_fout| _fout.close }
	[fout1, fout2, fout3].each{ |_fout| _fout.close }
end
desc "04.aln_extend_to_term"
task "04.aln_extend_to_term" do |t|
	dir  = ENV["dir"]
	Size = ENV["Size"].to_i
	odir = "#{dir}/result"

	Dir["#{odir}/03.aln.idt*.out"].each{ |fin|
		open("#{fin.split(".")[0..-2].join(".").gsub(/03\.aln/, "04.aln")}.extend.out", "w"){ |fout|
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
				fout.puts [a[0..1], idt, aln_len, mismatch, gap, q_start, q_stop, s_start, s_stop, a[10..-1]]*"\t"
			}
		}
	}
end
desc "05.make_circ_list_and_fasta"
task "05.make_circ_list_and_fasta" do |t|
	dir    = ENV["dir"]
	Idt    = ENV["Idt"].to_i
	Size   = ENV["Size"].to_i
	fin    = ENV["fin"]
	odir   = "#{dir}/result"
	fin0   = "#{odir}/04.aln.idt80.extend.out"
	fout   = open("#{odir}/05.circ.detected.out", "w")
	foutfa = open("#{odir}/05.circ.noTR.fasta", "w")   # ALL fasta (if circular, Terminal Redundancy is trimmed)

	node2remove_len = {}
	IO.readlines(fin0).select{ |l| l.chomp.split(/\t/)[2].to_f >= Idt }.each{ |l|
		fout.puts l
		a = l.chomp.split(/\t/)
		node = a[0]
		s_start, s_stop = a.values_at(8..9).map(&:to_i)
		node2remove_len[node] = s_stop - s_start + 1
	}
	IO.read(fin).split(/^>/)[1..-1].each{ |ent|
		node, *seq = ent.split("\n")
		seq = seq*""
		if remove_len = node2remove_len[node]
			seq = seq[0, seq.size - remove_len] # remove the last 'remove_len' seq
		end
		foutfa.puts [">"+node, seq.scan(/.{1,70}/)]
	}
end
