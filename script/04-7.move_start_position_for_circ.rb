
# fin0, fin1, fin, fout = ARGV
fin0, fin1, fgffs, fgff, fout = ARGV
fw = open(fout, "w")

def write_fasta(ent, newend, newlen, fw)
  if newend == newlen
    fw.puts ">"+ent
  else
    lab, seq = ent.split("\n")
    # p "newend:#{newend}"
    # p "size:#{seq[newend..-1].size}"
    fw.puts [">"+lab, seq[newend..-1] + seq[0..newend-1]]
  end
  return
end

def parse_and_write(ls, newlen, len, ent, fw)
  inter  = [] ## integergenic pos: [1..4, 200..210, ...]
  intra  = [] ## intragenic pos: [5..199, 211..300, ...]
  newend = 0 ## new end position to be determined (1-based)

  ls.each{ |l|
    id, _, _, start, stop = l.chomp.split("\t")
    start, stop = start.to_i, stop.to_i ## 1-based
    intra << (start..stop)
  }

  ## Is leftmost intergenic ?
  if (i = intra[0].begin) > 1
    inter << (1..i-1)
  end

  (2..intra.size).each{ |k|
    r1 = intra[k-2] ## intragenic range 1
    r2 = intra[k-1] ## intragenic range 2
    ## calc intergenic region between r1 and r2
    if r2.begin - r1.end >= 2 ## intergenic region exists
      inter << (r1.end+1..r2.begin-1)
    end
  }

  ## Is rightmost intergenic ?
  if (i = intra[-1].end) < len
    inter << (i+1..len)
  end

  inter.reverse.each{ |r| ### reverse
    if newlen < r.begin
      next
    elsif r.begin <= newlen and newlen <= r.end ## newlen is within integergenic region (no need to move)
      newend = newlen
      write_fasta(ent, newend, newlen, fw) 
      return
    elsif r.end < newlen
      newend = r.end
      write_fasta(ent, newend, newlen, fw) 
      return
    end
  }
end

### parse noTR fasta
id2ent = {}    ## noTR entry
id2newlen = {} ## noTR length
IO.read(fin0).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  id = lab.split(/\s+/)[0]
  id2ent[id] = ent
  id2newlen[id] = seq.join.gsub(/\s+/, "").size
}

### parse original fasta
id2len = {}    ## original length
IO.read(fin1).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  id = lab.split(/\s+/)[0]
  id2len[id] = seq.join.gsub(/\s+/, "").size
}

id2ls = Hash.new{ |h, i| h[i] = [] }
open(fgff, "w"){ |fw| ### write concat gff (sequence order is not the same as fasta)
  Dir[fgffs].each{ |fin|
    IO.readlines(fin).each{ |l|
      fw.puts l
      next if l =~ /^#/
      id = l.chomp.split("\t")[0]
      id2ls[id] << l
    }
  }
}

id2len.each_key{ |id|
  ls = id2ls[id]
  if ls.size > 0
    parse_and_write(ls, id2newlen[id], id2len[id], id2ent[id], fw)
  end
}

fw.close
