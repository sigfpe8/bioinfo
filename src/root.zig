//! By convention, root.zig is the root source file when making a library.
const std = @import("std");

const fas = @import("fasta.zig");
pub const Fasta = fas.Fasta;
pub const randomFastaArray    = fas.randomFastaArray;
pub const freeFastaArray      = fas.freeFastaArray;
pub const readFastaFile       = fas.readFastaFile;
pub const readFastaBuffer     = fas.readFastaBuffer;
pub const readFastaNoIdFile   = fas.readFastaNoIdFile;
pub const readFastaNoIdBuffer = fas.readFastaNoIdBuffer;
pub const writeFastaFile      = fas.writeFastaFile;
pub const fastaToSequences    = fas.fastaToSequences;
pub const getProtein          = fas.getProtein;

pub const readLine   = fas.readLine;
pub const readLines  = fas.readLines;
pub const freeLines  = fas.freeLines;
pub const printLines = fas.printLines;
pub const IntsReader = fas.IntsReader;

const utl = @import("utils.zig");
pub const binomial = utl.binomial;
pub const pascalsTriang = utl.pascalsTriang;
pub const PermIterator = utl.PermIterator;

const mot = @import("motif.zig");
pub const Motif = mot.Motif;

const assert = std.debug.assert;
const Io = std.Io;
const Reader = Io.Reader;
const Writer = Io.Writer;
const Allocator = std.mem.Allocator;
const Map = std.StringHashMap;

pub const BioError = error{
    DifferentLengths,
    TooManySequences,
};

/// Map one nucleotide to its complement
const dnaPairMap: [26]u8 = .{
    'T',    // A -> T
    'B',    // B
    'G',    // C -> G
    'D',    // D
    'E',    // E
    'F',    // F
    'C',    // G -> C
    'H',    // H
    'I',    // I
    'J',    // J
    'K',    // K
    'L',    // L
    'M',    // M
    'N',    // N
    'O',    // O
    'P',    // P
    'Q',    // Q
    'R',    // R
    'S',    // S
    'A',    // T -> A
    'U',    // U
    'V',    // V
    'W',    // W
    'X',    // X
    'Y',    // Y
    'Z',    // Z
};

/// Count of how many codons can encode one amino acid
pub const aacTable: [26]u8 = .{
//  3,    // .   UAA, UAG, UGA (Stop)
    4,    // A   GCU, GCC, GCA, GCG
    0,    // B
    2,    // C   UGU, UGC
    2,    // D   GAU, GAC
    2,    // E   GAA, GAG
    2,    // F   UUU, FFF
    4,    // G   GGU, GGC, GGA, GGG
    2,    // H   CAU, CAC
    3,    // I   AUU, AUC, AUA
    0,    // J
    2,    // K   AAA, AAG
    6,    // L   UUA, UUG, CUU, CUC, CUA, CUG
    1,    // M   AUG
    2,    // N   AAU, AAC
    0,    // O
    4,    // P   CCU, CCC, CCA, CGG
    2,    // Q   CAA, CAG
    6,    // R   CGU, CGC, CGA, CGG, AGA, AGG
    6,    // S   UCU, UCC, UCA, UCG, AGU, AGC
    4,    // T   ACU, ACC, ACA, ACG
    0,    // U
    4,    // V   GUU, GUC, GUA, GUG
    1,    // W   UGG
    0,    // X
    2,    // Y   UAU, UAC
    0,    // Z
};

pub const stopCodons: usize = 3;

/// Amino acid monoisotopic mass table
pub const aacMassTable: [26]f64 = .{
    71.03711 ,   // A
    0.0,         // B
    103.00919,   // C
    115.02694,   // D
    129.04259,   // E
    147.06841,   // F
    57.02146,    // G
    137.05891,   // H
    113.08406,   // I
    0.0,         // J
    128.09496,   // K
    113.08406,   // L
    131.04049,   // M
    114.04293,   // N
    0.0,         // O
    97.05276,    // P
    128.05858,   // Q
    156.10111,   // R
    87.03203,    // S
    101.04768,   // T
    0.0,         // U
    99.06841,    // V
    186.07931,   // W
    0.0,         // X
    163.06333,   // Y
    0.0,         // Z
};

const CodonMap = std.StaticStringMap(u8);

const rnaCodonMap = CodonMap.initComptime(.{
    .{"UUU", 'F'},    .{"CUU", 'L'},    .{"AUU", 'I'},    .{"GUU", 'V'},
    .{"UUC", 'F'},    .{"CUC", 'L'},    .{"AUC", 'I'},    .{"GUC", 'V'},
    .{"UUA", 'L'},    .{"CUA", 'L'},    .{"AUA", 'I'},    .{"GUA", 'V'},
    .{"UUG", 'L'},    .{"CUG", 'L'},    .{"AUG", 'M'},    .{"GUG", 'V'},
    .{"UCU", 'S'},    .{"CCU", 'P'},    .{"ACU", 'T'},    .{"GCU", 'A'},
    .{"UCC", 'S'},    .{"CCC", 'P'},    .{"ACC", 'T'},    .{"GCC", 'A'},
    .{"UCA", 'S'},    .{"CCA", 'P'},    .{"ACA", 'T'},    .{"GCA", 'A'},
    .{"UCG", 'S'},    .{"CCG", 'P'},    .{"ACG", 'T'},    .{"GCG", 'A'},
    .{"UAU", 'Y'},    .{"CAU", 'H'},    .{"AAU", 'N'},    .{"GAU", 'D'},
    .{"UAC", 'Y'},    .{"CAC", 'H'},    .{"AAC", 'N'},    .{"GAC", 'D'},
    .{"UAA", '.'},    .{"CAA", 'Q'},    .{"AAA", 'K'},    .{"GAA", 'E'},
    .{"UAG", '.'},    .{"CAG", 'Q'},    .{"AAG", 'K'},    .{"GAG", 'E'},
    .{"UGU", 'C'},    .{"CGU", 'R'},    .{"AGU", 'S'},    .{"GGU", 'G'},
    .{"UGC", 'C'},    .{"CGC", 'R'},    .{"AGC", 'S'},    .{"GGC", 'G'},
    .{"UGA", '.'},    .{"CGA", 'R'},    .{"AGA", 'R'},    .{"GGA", 'G'},
    .{"UGG", 'W'},    .{"CGG", 'R'},    .{"AGG", 'R'},    .{"GGG", 'G'}, 
});

/// Similar to codonMap{}, but includes the transciption T->U so that
/// the translation to amino acids can be done in a single step.
const dnaCodonMap = CodonMap.initComptime(.{
    .{"TTT", 'F'},    .{"CTT", 'L'},    .{"ATT", 'I'},    .{"GTT", 'V'},
    .{"TTC", 'F'},    .{"CTC", 'L'},    .{"ATC", 'I'},    .{"GTC", 'V'},
    .{"TTA", 'L'},    .{"CTA", 'L'},    .{"ATA", 'I'},    .{"GTA", 'V'},
    .{"TTG", 'L'},    .{"CTG", 'L'},    .{"ATG", 'M'},    .{"GTG", 'V'},
    .{"TCT", 'S'},    .{"CCT", 'P'},    .{"ACT", 'T'},    .{"GCT", 'A'},
    .{"TCC", 'S'},    .{"CCC", 'P'},    .{"ACC", 'T'},    .{"GCC", 'A'},
    .{"TCA", 'S'},    .{"CCA", 'P'},    .{"ACA", 'T'},    .{"GCA", 'A'},
    .{"TCG", 'S'},    .{"CCG", 'P'},    .{"ACG", 'T'},    .{"GCG", 'A'},
    .{"TAT", 'Y'},    .{"CAT", 'H'},    .{"AAT", 'N'},    .{"GAT", 'D'},
    .{"TAC", 'Y'},    .{"CAC", 'H'},    .{"AAC", 'N'},    .{"GAC", 'D'},
    .{"TAA", '.'},    .{"CAA", 'Q'},    .{"AAA", 'K'},    .{"GAA", 'E'},
    .{"TAG", '.'},    .{"CAG", 'Q'},    .{"AAG", 'K'},    .{"GAG", 'E'},
    .{"TGT", 'C'},    .{"CGT", 'R'},    .{"AGT", 'S'},    .{"GGT", 'G'},
    .{"TGC", 'C'},    .{"CGC", 'R'},    .{"AGC", 'S'},    .{"GGC", 'G'},
    .{"TGA", '.'},    .{"CGA", 'R'},    .{"AGA", 'R'},    .{"GGA", 'G'},
    .{"TGG", 'W'},    .{"CGG", 'R'},    .{"AGG", 'R'},    .{"GGG", 'G'}, 
});

const rnaCodonTable = [_][]const u8{
    "UUU", "F",      "CUU", "L",      "AUU", "I",      "GUU", "V",
    "UUC", "F",      "CUC", "L",      "AUC", "I",      "GUC", "V",
    "UUA", "L",      "CUA", "L",      "AUA", "I",      "GUA", "V",
    "UUG", "L",      "CUG", "L",      "AUG", "M",      "GUG", "V",
    "UCU", "S",      "CCU", "P",      "ACU", "T",      "GCU", "A",
    "UCC", "S",      "CCC", "P",      "ACC", "T",      "GCC", "A",
    "UCA", "S",      "CCA", "P",      "ACA", "T",      "GCA", "A",
    "UCG", "S",      "CCG", "P",      "ACG", "T",      "GCG", "A",
    "UAU", "Y",      "CAU", "H",      "AAU", "N",      "GAU", "D",
    "UAC", "Y",      "CAC", "H",      "AAC", "N",      "GAC", "D",
    "UAA", ".",      "CAA", "Q",      "AAA", "K",      "GAA", "E",
    "UAG", ".",      "CAG", "Q",      "AAG", "K",      "GAG", "E",
    "UGU", "C",      "CGU", "R",      "AGU", "S",      "GGU", "G",
    "UGC", "C",      "CGC", "R",      "AGC", "S",      "GGC", "G",
    "UGA", ".",      "CGA", "R",      "AGA", "R",      "GGA", "G",
    "UGG", "W",      "CGG", "R",      "AGG", "R",      "GGG", "G", 
};

pub fn makeRNACodonMap(allocator: Allocator) !Map(u8) {
    var map = Map(u8).init(allocator);
    var i: usize = 0;

    while (i < rnaCodonTable.len - 1) : (i += 2) {
        const rna = rnaCodonTable[i];
        const aac = rnaCodonTable[i+1][0];
        try map.put(rna,aac);
    }

    return map;
}

/// For each amino acid, determine how many codons may encode it
pub fn makeAminoAcidTable() void {
    var table = [_]u8{0} ** 26;
    var i: usize = 0;

    while (i < rnaCodonTable.len - 1) : (i += 2) {
        // const rna = rnaCodonTable[i];
        const aac = rnaCodonTable[i+1][0];
        if (aac >= 'A' and aac <= 'Z') {
            table[aac - 'A'] += 1;
        }
    }

    std.debug.print("pub const aacTable: [26]u8 = .{{\n", .{});
    i = 0;
    while (i < table.len) : (i += 1) {
        const c: u8 = @truncate(i + 'A');
        std.debug.print("    {d},    // {c}\n", .{table[i], c});
    }
    std.debug.print("}};\n", .{});
}

/// Check if the sequence is valid (only contains A, C, G, T).
pub fn validSeq(seq: []const u8) bool {
    for (seq) |c| {
        if (c != 'A' and c != 'C' and c != 'G' and c != 'T') {
            return false;
        }
    }
    return true;
}

/// Generate a random DNA sequence of the given length.
/// Caller owns the returned string.
pub fn randomSeq(allocator: Allocator, len: usize) ![]u8 {
    const seq = try allocator.alloc(u8, len);
    var prng = std.Random.DefaultPrng.init(0); 
    const rand = prng.random();
    for (seq) |*c| {
        const r = rand.uintLessThan(usize, 4);
        switch (r) {
            0 => c.* = 'A',
            1 => c.* = 'C',
            2 => c.* = 'G',
            3 => c.* = 'T',
            else => unreachable,
        }
    }
    return seq;
}

/// Count the number of each base in the sequence and return a tuple of (A, C, G, T).
pub fn countBases(seq: []const u8) struct { usize, usize, usize, usize } {
    var nA: usize = 0;
    var nC: usize = 0;
    var nG: usize = 0;
    var nT: usize = 0;

    for (seq) |c| {
        switch (c) {
            'A' => nA += 1,
            'C' => nC += 1,
            'G' => nG += 1,
            'T' => nT += 1,
            else => {}, // Might be \n
        }
    }

    return .{ nA, nC, nG, nT };
}

/// Return the percentage of GC content in the sequence.
pub fn gcContent(seq: []const u8) f64 {
    const nA, const nC, const nG, const nT = countBases(seq);
    const total = nA + nC + nG + nT;
    if (total == 0) {
        return 0.0;
    }
    return (@as(f64,@floatFromInt(nC + nG)) / @as(f64, @floatFromInt(total))) * 100.0;
}

/// Return the complement of the given DNA sequence.
/// Caller owns the returned string.
pub fn complement(allocator: Allocator, seq: []const u8) ![]u8 {
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        comp[i] = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
    }
    return comp;
}

/// Return the reverse complement of the given DNA sequence.
/// Caller owns the returned string.
pub fn revComplement(allocator: Allocator, seq: []const u8) ![]u8 {
    const len = seq.len - 1;
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        comp[len - i] = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
    }
    return comp;
}

/// Check if a sequence is a "Reverse Palindrome", that is, if
/// its reverse complement is equal to itself.
pub fn isRevPalindrome(seq: []const u8) bool {
    const len = seq.len;
    if (len == 0 or len % 2 == 1) {
        return false;
    }

    var i: usize = 0;
    while (i < len / 2) : (i += 1) {
        const c = seq[i];
        const r = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
        if (seq[len - i - 1] != r) {
            return false;
        }
    }

    return true;
}

/// Transcribe the given DNA sequence to RNA (replace T with U).
/// The caller is responsible for freeing the memory of the returned slice.
/// Caller owns the returned string.
pub fn transcribe(allocator: Allocator, seq: []const u8) ![]u8 {
    var rna = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        rna[i] = if (c == 'T') 'U' else c;
    }
    return rna;
}

/// Translate an RNA sequence into an amino acid sequence (a protein).
/// If the last codon encodes a stop, remove it from the output.
/// Caller owns the returned string.
pub fn translateRNA(allocator: Allocator, rna: []const u8) ![]u8 {
    assert(rna.len % 3 == 0);
    var len = rna.len;

    // If the input sequence ends in a 'stop' codon, ignore it
    if (rnaCodonMap.get(rna[len-3..]) == '.') {
        len -= 3;
    }

    var prot = try allocator.alloc(u8, len / 3);

    var i: usize = 0;
    while (i < len) : (i += 3) {
        const cod = rna[i..i+3];             // The 3-letter codon
        const aac = rnaCodonMap.get(cod).?;  // Corresponding amino acid
        prot[i / 3] = aac;
    }

    return prot;
}

/// Translate an DNA sequence into an amino acid sequence (a protein).
/// Similar to translateRNA() but goes directly from DNA to protein.
/// Caller owns the returned string.
pub fn translateDNA(allocator: Allocator, dna: []const u8) ![]u8 {
    assert(dna.len % 3 == 0);
    var len = dna.len;

    // If the input sequence ends in a 'stop' codon, ignore it
    if (dnaCodonMap.get(dna[len-3..]) == '.') {
        len -= 3;
    }

    var prot = try allocator.alloc(u8, len / 3);

    var i: usize = 0;
    while (i < len) : (i += 3) {
        const cod = dna[i..i+3];             // The 3-letter codon
        const aac = dnaCodonMap.get(cod).?;  // Corresponding amino acid
        prot[i / 3] = aac;
    }

    return prot;
}

/// Given a DNA sequence, remove all the introns in it and then translate it
/// to the encoded sequence of amino acids (a protein).
/// Caller owns the returned string.
pub fn spliceDNA(allocator: Allocator, seq: []const u8, introns: [][]u8) ![]u8 {
    var buffer = try allocator.dupe(u8, seq);   // Work buffer
    defer allocator.free(buffer);
    var dna = buffer;

    for (introns) |intron| {
        const len = cutIntron(dna, intron);
        dna = buffer[0..len];
    }

    return translateDNA(allocator, dna);
}

/// Cut off all occurrences of the intron in this sequence.
/// Return the length of the shortened input sequence.
/// Cuts are made in place.
fn cutIntron(buf: []u8, intron: []u8) usize {
    var base: usize = 0;
    var len: usize = buf.len;

    while (std.mem.find(u8, buf[base..], intron)) |pos| {
        const end = pos + intron.len;                // Where this intron ends
        const lft = len - end;                       // How many nucleotides left after the intron
        @memmove(buf[pos..pos+lft], buf[end..len]);  // Move everything after the intron
        len -= intron.len;
        base = pos;
    }

    return len;
}


/// Return the Hamming distance between two sequences, that is,
/// the number of corresponding nucleotides that are different.
/// If the sequences have different lengths, return -1.
pub fn hammingDist(seq1: []const u8, seq2: []const u8) isize {
    if (seq1.len != seq2.len) {
        return -1;
    }

    var hamd: isize = 0;

    for (0..seq1.len) |i| {
        if (seq1[i] != seq2[i]) {
            hamd += 1;
        }
    }

    return hamd;
}

pub fn profile(allocator: Allocator, seqs: [][]u8) ![][]u8 {
    // Since we use u8's to count the # of each nucleotide in a certain
    // position, the maximum # of sequences we can compare is 255. This
    // seems reasonable for the current tests buf if a larger number is
    // required, we can make the type a comptime parameter and use bigger
    // integers when necessary.
    if (seqs.len > 255) {
        return BioError.TooManySequences;
    }

    // Make sure all the sequences have the same length
    const len = seqs[0].len;
    for (seqs[1..]) |seq| {
        if (seq.len != len) {
            return BioError.DifferentLengths;
        }
    }

    // Allocate the profile matrix and initialize it with 0's.
    // Each line has the same length as the incoming DNA sequences.
    //    A:  0 0 0 ... 0
    //    C:  0 0 0 ... 0
    //    G:  0 0 0 ... 0
    //    T:  0 0 0 ... 0
    var prof = try allocator.alloc([]u8,4);
    for (0..4) |i| {
        const cnt = try allocator.alloc(u8, len);
        @memset(cnt, 0);
        prof[i] = cnt;
    }

    // Scan all the incoming DNA sequences
    for (seqs) |seq| {
        // Scan all positions
        for (seq, 0..) |c, p| {
            const n: usize = switch (c) {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                else => unreachable,
            };
            // Increment counter for corresponding nucleotide and position
            prof[n][p] += 1;
        }
    }

    return prof;
}

pub fn freeProfile(allocator: Allocator, prof: [][]u8) void {
    for (prof) |cnts| {
        allocator.free(cnts);
    }
    allocator.free(prof);
}

pub fn consensus(allocator: Allocator, prof: [][]u8) ![]u8 {
    const len = prof[0].len;
    var cons = try allocator.alloc(u8, len);

    for (0..len) |i| {
        var max = prof[0][i];
        var idx: usize = 0;
        for (1..4) |j| {
            if (prof[j][i] > max) {
                max = prof[j][i];
                idx = j;
            }
        }
        cons[i] = switch (idx) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            else => unreachable,
        };
    }

    return cons;
}

/// Return the total mass (in Da) of a protein string
pub fn proteinMass(pro: []const u8) f64 {
    var mass: f64 = 0.0;

    for (pro) |a| {
        if (a >= 'A' and a <= 'Z') {
            mass += aacMassTable[a - 'A'];
        }
    }

    return mass;
}

test "validSeq" {
    try std.testing.expect(validSeq("ACGT"));
    try std.testing.expect(!validSeq("ACGTX"));
}

test "randomSeq" {
    const allocator = std.testing.allocator;
    const seq = try randomSeq(allocator, 100);
    defer allocator.free(seq);
    try std.testing.expect(validSeq(seq));
}

test "countBases" {
    const seq = "AAGCTAGC";
    const nA, const nC, const nG, const nT = countBases(seq);
    try std.testing.expectEqual(3, nA); // A
    try std.testing.expectEqual(2, nC); // C
    try std.testing.expectEqual(2, nG); // G
    try std.testing.expectEqual(1, nT); // T
}

test "gcContent" {
    const seq = "AAGCTAGC";
    const gc = gcContent(seq);
    try std.testing.expectEqual(50.0, gc); // 4 out of 8 are G or C
}

test "complement" {
    const allocator = std.testing.allocator;
    const seq = "ACGT\nTGCA";
    const comp = try complement(allocator, seq);
    defer allocator.free(comp);
    try std.testing.expectEqualSlices(u8, "TGCA\nACGT", comp);
}

test "revComplement" {
    const allocator = std.testing.allocator;
    const seq = "ACGT\nTGCA";
    const comp = try revComplement(allocator, seq);
    defer allocator.free(comp);
    try std.testing.expectEqualSlices(u8, "TGCA\nACGT", comp);
}

test "transcribe" {
    const allocator = std.testing.allocator;
    const seq = "ACGT";
    const rna = try transcribe(allocator, seq);
    defer allocator.free(rna);
    try std.testing.expectEqualSlices(u8, "ACGU", rna);
}

test "hammingDist" {
    try std.testing.expectEqual(hammingDist("CAAACGTT","GAAACGT"), -1);
    try std.testing.expectEqual(hammingDist("CAAACGTT","CAAACGTT"), 0);
    try std.testing.expectEqual(hammingDist("CAAACGTT","CAAACGTC"), 1);
    try std.testing.expectEqual(hammingDist("CAAACGTT","GAAACGTC"), 2);
}

test "Motif" {
    const gpa = std.testing.allocator;
    var motif = Motif.init(gpa);
    defer motif.deinit();
    // try motif.encode("ABC[DEF]G{HIJKL}");

    try motif.encode("ABC");
    try std.testing.expect(motif.match("ABC"));     // Exact match
    try std.testing.expect(motif.match("ABCD"));    // Pattern shorter than text
    try std.testing.expect(!motif.match("AB"));     // Text shorter than pattern (fail)

    try motif.encode("A[BCD]E");
    try std.testing.expect(motif.match("ABE"));
    try std.testing.expect(motif.match("ACE"));
    try std.testing.expect(motif.match("ADE"));
    try std.testing.expect(!motif.match("AFE"));

    try motif.encode("A{BCD}E");
    try std.testing.expect(!motif.match("ABE"));
    try std.testing.expect(!motif.match("ACE"));
    try std.testing.expect(!motif.match("ADE"));
    try std.testing.expect(motif.match("AFE"));

    try motif.encode("ABC[DEF]G{HIJKL}");
    try std.testing.expect(motif.match("ABCEGM"));
    try std.testing.expect(motif.match("ABCFGYZ"));
    try std.testing.expect(!motif.match("ABCFYZ"));
}