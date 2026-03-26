//! By convention, root.zig is the root source file when making a library.
const std = @import("std");
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
/// The caller is responsible for freeing the memory.
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
pub fn complement(allocator: Allocator, seq: []const u8) ![]u8 {
    var comp = try allocator.alloc(u8, seq.len);
    for (seq, 0..) |c, i| {
        var pair = c;
        if (c >= 'A' and c <= 'Z') {
            pair = dnaPairMap[c - 'A'];
        }
        comp[i] = pair;
    }
    return comp;
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