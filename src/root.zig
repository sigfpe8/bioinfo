//! By convention, root.zig is the root source file when making a library.
const std = @import("std");
const assert = std.debug.assert;
const Io = std.Io;
const Reader = Io.Reader;
const Writer = Io.Writer;
const Allocator = std.mem.Allocator;

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
/// The caller is responsible for freeing the memory of the returned slice.
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
pub fn complement(allocator: Allocator, seq: []const u8) ![]u8 {
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        comp[i] = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
    }
    return comp;
}

/// Return the reverse complement of the given DNA sequence.
/// The caller is responsible for freeing the memory of the returned slice.
pub fn revComplement(allocator: Allocator, seq: []const u8) ![]u8 {
    const len = seq.len - 1;
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        comp[len - i] = if (c >= 'A' and c <= 'Z') dnaPairMap[c - 'A'] else c;
    }
    return comp;
}

/// Transcribe the given DNA sequence to RNA (replace T with U).
/// The caller is responsible for freeing the memory of the returned slice.
pub fn transcribe(allocator: Allocator, seq: []const u8) ![]u8 {
    var rna = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        rna[i] = if (c == 'T') 'U' else c;
    }
    return rna;
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


// Support for FASTA files

// This struct represents a single sequence in a FASTA file.
// A file can contain multiple sequences.
const Fasta = struct {
    seqID: []const u8,
    seq: []const u8,
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

pub fn readFastaFile(io: Io, allocator: Allocator, fname: []const u8) ![]Fasta {
    const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
    defer file.close(io);
    var buffer: [1024]u8 = undefined;
    var file_reader: Io.File.Reader = .init(file, io, &buffer);

    return try readFasta(allocator, &file_reader.interface);
}

pub fn readFastaBuffer(allocator: Allocator, buffer: []const u8) ![]Fasta {
    var reader: Reader = .fixed(buffer);

    return try readFasta(allocator, &reader);
}

pub fn readFasta(allocator: Allocator, r: *Reader) ![]Fasta {
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

const maxFastaSeqLine = 80;

pub fn writeFastaFile(io: Io, seqs: []const Fasta, fname: []const u8) !void {
    const file = try Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);
    var buffer: [1024]u8 = undefined;
    var file_writer: Io.File.Writer = .init(file, io, &buffer);

    return try writeFasta(seqs, &file_writer.interface);
}

pub fn writeFasta(seqs: []const Fasta, w: *Writer) !void {
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