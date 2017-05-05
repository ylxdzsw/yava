include("OhMyJulia.jl")
include("BioDataStructures.jl")
include("Falcon.jl")
include("Fire.jl")

using OhMyJulia
using BioDataStructures
using Falcon
using Fire

function load_homo()
    result = Dict{String, String}()
    for line in eachline("homo.data")
        s = findlast(line, '\t')
        result[line[1:s-1]] = line[s+1:end-1]
    end
    result
end

function load_freq()
    f90, f98 = fill(">20%", 1000, 240), fill(">20%", 1000, 240)
    for (d, line) in enumerate(eachline(split, "freq.data"))
        a = parse(Int, car(line)[1:end-1])
        for i in 1:a
            f90[d, i] = line[2i]
            f98[d, i] = line[2i+1]
        end
    end
    f90, f98
end

@main function anno()
    homo = load_homo()
    f90, f98 = load_freq()

    nall, nhomo, na0 = 0, 0, 0

    STDOUT << chomp(readline(STDIN)) << " | yavainfo: homoseq, 90% confidence freq, 98% confidence freq\n"
    for line in eachline(chomp, STDIN)
        STDOUT << line
        chr, pos, ref, alt, info = split(line, '\t')[[1, 2, 4, 5, end]]
        nall += 1

        STDOUT << '\t'
        var = join((chr, pos, ref, alt), '\t')
        if var in keys(homo)
            nhomo += 1
            STDOUT << homo[var]
        else
            STDOUT << 0
        end

        mrbam = map(x->parse(Int, x), split(info[findlast(info, ':')+1:end], ','))
        d, a = sum(mrbam), sum(mrbam[7:end])

        if a == 0
            na0 += 1
        end

        while d > 1000
            d >>= 1
            a >>= 1
        end

        STDOUT << '\t'
        if a == 0 || d == 0
            STDOUT << "0%\t0%"
        elseif a > 240
            STDOUT << ">20%\t>20%"
        else
            STDOUT << f90[d, a] << '\t' << f98[d, a]
        end

        STDOUT << '\n'
    end

    STDERR << """
    ==== yava ====
    yava finished succesfully
    total number of variants annotated: $nall
    variants that can be caused by similar sequence: $nhomo
    variants that has no read support: $na0
    ==== yava ====
    """ 
end