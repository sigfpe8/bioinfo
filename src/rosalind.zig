//! Rosalind problems and solutions
//! Bioinformatics stronghold: https://rosalind.info/problems/list-view/

const std = @import("std");
const bio = @import("root.zig");
const Allocator = std.mem.Allocator;
const print = std.debug.print;

var gpa: Allocator = undefined;
var io: std.Io = undefined;

pub fn main(init: std.process.Init) !void {
    gpa = init.gpa;
    io = init.io;
    // const seq = "GATTCGTGGACAGTGCATAAACCGTCGCGTTGTCTGTTCAGCTCGGGGCAGGCGATCCCATAGTACCTCAGGAGCTCGGACGGATGTCACCCAGGTTAATGAAGAGAGTCTAACACCTATTATGCAGGGCGATGGCCCATAATGGGTTCGCTTTACCGCCTGCGGATGCATGGGGTAGCGAGTAGAAAGTTCATTGGGTCCCAGGAGACCTACGACTGAATGCAGACGCGATTTTCCGCAGCTAACATCTGGACGTTAATTCCGCCGGGGCGTTGTCAATAGTGCACTCAACCCCGCTGAAACACGAACGCAACGCTGGTCTGCAACTTAAAAGAGCGCTCTTAACGGGCTTACTGTATCCCGGGATGCAATTGTCCCGTGTAGTCCCGCCCTCCTTAGTTATATATTGCAGCAATCAAGGGCGCTAGTTCACTATGCTTGGCACGTTATAAACTTAAGTTTTACTGTTCGAGTTCGGTCGTGCGGGGAAAGCGTGGCGTGACTCATCTGTTTAATGGGGTCCAGCTATAATCCAGTAATCTAATAAGTCAGCCGGTCCTTGGTCTCTTCCACCCTCCTTGTTTTTCGTGCGCTAGTAACACATCTATTCAGGCAGCCGACCAGGTCCCTCCTCCCCCAGTTCATTGCCGACACTCGTCACTGCGCGGCGCCTGGTCCACGCTACATCGATGTAATTCATGCACTGTGAGGTGTACCTCGCACATGAATCTAACCACTGCGAAAGTGGGTGGTGCGCAATACGTGTGCGGCTCGCCAGGTAGGTACCATGGATGACTTCCTAGTGTGGAGCTGCCTCAACAATCATGTTAAAG";
    // problemDNA(seq);
    // try problemRNA(seq);
    //try problemREVC(seq);
    // try problemGC();
    // problemFIB(28, 5);
    problemHAMM(
        "GAGCCTACTAACGGGAT",
        "CATCGTAATGACGGCCT"
    );
}

/// DNA - Counting DNA Nucleotides
pub fn problemDNA(seq: []const u8) void {
    const nA, const nC, const nG, const nT = bio.countBases(seq);
    print("{d} {d} {d} {d}\n", .{nA,nC,nG,nT});
}

/// RNA - Transcribing DNA into RNA
pub fn problemRNA(seq: []const u8) !void {
    const rna = try bio.transcribe(gpa, seq);
    print("{s}\n", .{rna});
    gpa.free(rna);
}

/// REVC - Complementing a Strand of DNA
pub fn problemREVC(seq: []const u8) !void {
    const revC = try bio.revComplement(gpa, seq);
    print("{s}\n", .{revC});
    gpa.free(revC);
}

/// GC - Computing GC content
pub fn problemGC() !void {
    // const fname = "~/Downloads/rosalind_gc.txt";
    const fname = "/Users/cordeiro/Downloads/rosalind_gc.txt";
    const seqs = try bio.readFastaFile(io, gpa, fname);

    if (seqs.len > 0) {
        var gc = bio.gcContent(seqs[0].seq);
        var id: usize = 0;

        for (seqs[1..], 1..) |s, i| {
            const x = bio.gcContent(s.seq);
            if (x > gc) {
                gc = x;
                id = i;
            }
            // print("{d}: {s}\n", .{ i, s.seqID });
            // print("{s}\n\n", .{s.seq});
        }
        print("{s}\n", .{seqs[id].seqID});
        print("{d:.6}\n", .{gc});
    }

    bio.freeFastaArray(gpa, seqs);
}

pub fn problemFIB(n: usize, k: usize) void {
    var fn2: usize = 1;
    var fn1: usize = 1;
    var f: usize = 1;
    var i: usize = 3;

    while (i <= n) : (i += 1) {
        // Recurrence relation: Fn = Fn-1 + k * Fn-2
        f = fn1 + k * fn2;
        fn2 = fn1;
        fn1 = f;
    }

    print("{d}\n", .{f});
}

pub fn problemHAMM(seq1: []const u8, seq2: []const u8) void {
    const hamd = bio.hammingDist(seq1, seq2);

    print("{d}\n", .{hamd});
}