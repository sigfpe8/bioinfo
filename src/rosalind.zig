//! Rosalind problems and solutions
//! Bioinformatics stronghold: https://rosalind.info/problems/list-view/

const std = @import("std");
const bio = @import("root.zig");
const Allocator = std.mem.Allocator;
const print = std.debug.print;

var gpa: Allocator = undefined;

pub fn main(init: std.process.Init) !void {
    gpa = init.gpa;
    const seq = "GATTCGTGGACAGTGCATAAACCGTCGCGTTGTCTGTTCAGCTCGGGGCAGGCGATCCCATAGTACCTCAGGAGCTCGGACGGATGTCACCCAGGTTAATGAAGAGAGTCTAACACCTATTATGCAGGGCGATGGCCCATAATGGGTTCGCTTTACCGCCTGCGGATGCATGGGGTAGCGAGTAGAAAGTTCATTGGGTCCCAGGAGACCTACGACTGAATGCAGACGCGATTTTCCGCAGCTAACATCTGGACGTTAATTCCGCCGGGGCGTTGTCAATAGTGCACTCAACCCCGCTGAAACACGAACGCAACGCTGGTCTGCAACTTAAAAGAGCGCTCTTAACGGGCTTACTGTATCCCGGGATGCAATTGTCCCGTGTAGTCCCGCCCTCCTTAGTTATATATTGCAGCAATCAAGGGCGCTAGTTCACTATGCTTGGCACGTTATAAACTTAAGTTTTACTGTTCGAGTTCGGTCGTGCGGGGAAAGCGTGGCGTGACTCATCTGTTTAATGGGGTCCAGCTATAATCCAGTAATCTAATAAGTCAGCCGGTCCTTGGTCTCTTCCACCCTCCTTGTTTTTCGTGCGCTAGTAACACATCTATTCAGGCAGCCGACCAGGTCCCTCCTCCCCCAGTTCATTGCCGACACTCGTCACTGCGCGGCGCCTGGTCCACGCTACATCGATGTAATTCATGCACTGTGAGGTGTACCTCGCACATGAATCTAACCACTGCGAAAGTGGGTGGTGCGCAATACGTGTGCGGCTCGCCAGGTAGGTACCATGGATGACTTCCTAGTGTGGAGCTGCCTCAACAATCATGTTAAAG";
    // problemDNA(seq);
    // try problemRNA(seq);
    try problemREVC(seq);
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