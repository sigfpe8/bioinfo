//! By convention, root.zig is the root source file when making a library.
const std = @import("std");

const utl = @import("utils.zig");
pub const binomial = utl.binomial;
pub const readLine = utl.readLine;
pub const readLines = utl.readLines;
pub const freeLines = utl.freeLines;
pub const printLines = utl.printLines;
pub const IntsReader = utl.IntsReader;
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
    FailedToLoadProtein,
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
    3,    // I   AUU, AUC
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

const codonMap = CodonMap.initComptime(.{
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
/// The caller is responsible for freeing the memory of the returned slice.
/// Caller owns the returned string.
pub fn complement(allocator: Allocator, seq: []const u8) ![]u8 {
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        comp[i] = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
    }
    return comp;
}

/// Return the reverse complement of the given DNA sequence.
/// The caller is responsible for freeing the memory of the returned slice.
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
pub fn translate(allocator: Allocator, rna: []const u8) ![]u8 {
    assert(rna.len % 3 == 0);
    var len = rna.len;

    // If the input sequence ends in a 'stop' codon, ignore it
    if (codonMap.get(rna[len-3..]) == '.') {
        len -= 3;
    }

    var prot = try allocator.alloc(u8, len / 3);

    var i: usize = 0;
    while (i < len) : (i += 3) {
        const cod = rna[i..i+3];            // The 3-letter codon
        const aac = codonMap.get(cod).?;    // Corresponding amino acid
        prot[i / 3] = aac;
    }

    return prot;
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

// Support for FASTA files

// This struct represents a single sequence in a FASTA file.
// A file can contain multiple sequences.
pub const Fasta = struct {
    seqID: []u8,
    seq: []u8,
};

pub fn randomFastaArray(allocator: Allocator, n: usize) ![]Fasta {
    var prng = std.Random.DefaultPrng.init(0); 
    const rand = prng.random();

    var list: std.ArrayList(Fasta) = .empty;
    var i: usize = 1;

    while (i <= n) : (i += 1) {
        const seqId = try std.fmt.allocPrint(allocator, "Sequence {d}", .{i});
        const len = rand.uintLessThan(usize, 100) + 50;
        const seq = try randomSeq(allocator, len);
        try list.append(allocator, .{ .seqID = seqId, .seq = seq });
    }

    return list.toOwnedSlice(allocator);
}

pub fn freeFastaArray(allocator: Allocator, array: []Fasta) void {
    for (array) |s| {
        allocator.free(s.seqID);
        allocator.free(s.seq);
    }
    allocator.free(array);
}

/// Reads a file in FASTA format and returns an array of Fasta structures.
/// Joins all the lines from the same sequence into a single string.
/// Caller must call freeFastaArray() to dispose of the array.
pub fn readFastaFile(io: Io, allocator: Allocator, fname: []const u8) ![]Fasta {
    const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
    defer file.close(io);
    var buffer: [1024]u8 = undefined;
    var file_reader: Io.File.Reader = .init(file, io, &buffer);

    return try readFasta(allocator, &file_reader.interface);
}

/// Reads a buffer in FASTA format and returns an array of Fasta structures.
/// Joins all the lines from the same sequence into a single string.
/// Caller must call freeFastaArray() to dispose of the array.
pub fn readFastaBuffer(allocator: Allocator, buffer: []const u8) ![]Fasta {
    var reader: Reader = .fixed(buffer);

    return try readFasta(allocator, &reader);
}

/// Main code for reading FASTA sequences
fn readFasta(allocator: Allocator, r: *Reader) ![]Fasta {
    var list: std.ArrayList(Fasta) = .empty;
    var seqId: ?[]u8 = null;
    var seq: std.ArrayList(u8) = .empty;

    while (true) {
        const line = r.takeDelimiterInclusive('\n') catch |err| {
            switch (err) {
                error.EndOfStream => break,
                else => return err,
            }
        };
        if (line.len == 1) {    // Ignore empty lines (just \n) 
            continue;
        }
        if (line[0] == '>') {   // New sequence ID
            if (seqId) |id| { // If there's a previous sequence, store it
                // Store previous sequence
                try list.append(allocator, .{
                            .seqID = id,
                            .seq = try seq.toOwnedSlice(allocator) });
            }
            // Mark new sequence ID
            seqId = try allocator.dupe(u8, line[1..line.len - 1]);
            seq = .empty;
        } else {
            // Concatenate lines in the same sequence
            try seq.appendSlice(allocator, line[0..line.len - 1]); // Trim \n
        }
    }

    // Include last sequence
    if (seqId) |id| {
        try list.append(allocator, .{
                    .seqID = id,
                    .seq = try seq.toOwnedSlice(allocator) });
    }

    return list.toOwnedSlice(allocator);
}

/// Similar to readFastaFile(), but ignores the IDs.
/// Returns an array of lines like readLines().
/// Caller must call freeLines() to dispose of the array.
pub fn readFastaNoIdFile(io: Io, allocator: Allocator, fname: []const u8) ![][]u8 {
    const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
    defer file.close(io);
    var buffer: [1024]u8 = undefined;
    var file_reader: Io.File.Reader = .init(file, io, &buffer);

    return try readFastaNoId(allocator, &file_reader.interface);
}

/// Similar to readFastaBuffer(), but ignores the IDs.
/// Returns an array of lines like readLines().
/// Caller must call freeLines() to dispose of the array.
pub fn readFastaNoIdBuffer(allocator: Allocator, buffer: []const u8) ![][]u8 {
    var reader: Reader = .fixed(buffer);

    return try readFastaNoId(allocator, &reader);
}

/// Similar to readFasta(), but ignores the IDs.
/// Returns an array of lines like readLines().
/// Caller must call freeLines() to dispose of the array.
fn readFastaNoId(allocator: Allocator, r: *Reader) ![][]u8 {
    var lines: std.ArrayList([]u8) = .empty;
    var seq: std.ArrayList(u8) = .empty;
    var pend: bool = false;     // True if there's a sequence pending
    var line: []u8 = undefined;

    while (true) {
        line = r.takeDelimiterInclusive('\n') catch |err| {
            switch (err) {
                error.EndOfStream => break,
                else => return err,
            }
        };
        if (line.len == 1) {    // Ignore empty lines (just \n) 
            continue;
        }
        if (line[0] == '>') {   // New sequence ID
            if (pend) { // If there's a previous sequence, store it
                // Store previous sequence
                line = try seq.toOwnedSlice(allocator);
                try lines.append(allocator, line);
            }
            // Mark new sequence ID
            pend = true;
            seq = .empty;
        } else {
            // Concatenate lines in the same sequence
            try seq.appendSlice(allocator, line[0..line.len - 1]); // Trim \n
        }
    }

    // Include last sequence
    if (pend) {
        line = try seq.toOwnedSlice(allocator);
        try lines.append(allocator, line);
    }

    return lines.toOwnedSlice(allocator);
}

const maxFastaSeqLine = 80;

pub fn writeFastaFile(io: Io, seqs: []const Fasta, fname: []const u8) !void {
    const file = try Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);
    var buffer: [1024]u8 = undefined;
    var file_writer: Io.File.Writer = .init(file, io, &buffer);

    return try writeFasta(seqs, &file_writer.interface);
}

fn writeFasta(seqs: []const Fasta, w: *Writer) !void {
    // >SeqID
    // seq
    for (seqs) |s| {
        try w.print(">{s}\n", .{s.seqID});
        var p = s.seq;
        // Print sequence in chunks of maxFastaSeqLine characters
        while (p.len > 0) {
            const l = @min(p.len, maxFastaSeqLine);
            try w.print("{s}\n", .{p[0..l]});
            p = p[l..]; // Next chunk
        }
    }
    try w.flush();
}

/// Duplicate an array of Fasta sequences as a simple array of sequences
pub fn fastaToSequences(allocator: Allocator, seqs: []Fasta) ![][]u8 {
    var lines = try allocator.alloc([]u8, seqs.len);

    for (seqs, 0..) |seq, i| {
        lines[i] = try allocator.dupe(u8, seq.seq);
    }

    return lines;
}

/// Get a protein fasta file
/// First try to get it from the cache (datasets/uniprot).
/// If not already there, fetch it from the site uniprot.org and save it in the cache.
pub fn getProtein(io: Io, gpa: Allocator, prot_name: []const u8) ![]Fasta {
    const fname = try std.fmt.allocPrint(gpa, "datasets/uniprot/{s}.fasta", .{prot_name});
    defer gpa.free(fname);

    var fasta: []Fasta = undefined;

    fasta = readFastaFile(io, gpa, fname) catch not_cached: {
        // Looks like this protein is not cached; fetch it from the web
        const prot = fetchProtein(io, gpa, prot_name) catch not_fetched: {
            // Fetch from the Web failed
            // Try simplified name, say P22457 instead of P22457_FA7_BOVIN
            if (std.mem.findScalar(u8, prot_name, '_')) |pos| {
                const pro = try fetchProtein(io, gpa, prot_name[0..pos]);
                break :not_fetched pro;
            } else {
                return BioError.FailedToLoadProtein;
            }
        };
        defer gpa.free(prot);
        // Save the file contents for future use
        try cacheProtein(io, gpa, prot_name, prot);
        const fetched = try readFastaBuffer(gpa, prot);
        break :not_cached fetched;
    };

    return fasta;
}

/// Fetch a protein .fasta file from uniprot.org
/// Return the contents of the file as a string so that it can be
/// easily cached and then be used to generate the []Fasta array.
pub fn fetchProtein(io: Io, gpa: Allocator, prot_name: []const u8) ![]u8 {
    var client: std.http.Client = .{
        .allocator = gpa,
        .io = io,
    };
    defer client.deinit();

    var result_body = std.Io.Writer.Allocating.init(gpa);
    defer result_body.deinit();

    const prot_site = "https://www.uniprot.org/uniprot/";
    const url = try std.fmt.allocPrint(gpa,"{s}{s}.fasta", .{prot_site,prot_name});
    defer gpa.free(url);

    const response = try client.fetch(.{
        .location = .{ .url = url },
        .response_writer = &result_body.writer,
    });

    if (response.status.class() != .success) {
        return BioError.FailedToLoadProtein;
    }

    return result_body.toOwnedSlice();
}

/// Save a just fetched protein in the cache directory.
/// `prot` is the full contents of the fasta file as a string.
fn cacheProtein(io: Io, gpa: Allocator, prot_name: []const u8, prot: []const u8) !void {
    const fname = try std.fmt.allocPrint(gpa, "datasets/uniprot/{s}.fasta", .{prot_name});
    defer gpa.free(fname);

    const file = try Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [1024]u8 = undefined;
    var file_writer: Io.File.Writer = .init(file, io, &buffer);
    var writer = &file_writer.interface;

    try writer.print("{s}\n", .{prot});
    try writer.flush();
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