const std = @import("std");
const bio = @import("bioinfo");
const Allocator = std.mem.Allocator;
const print = std.debug.print;

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;
    const seq = try bio.randomSeq(gpa, 50);
    defer gpa.free(seq);
    std.debug.print("Random sequence: {s}\n", .{seq});

    const fname = "example.fas";
    const seqs = try bio.readFastaFile(io, gpa, fname);

    print("Sequences in file '{s}':\n", .{fname});
    for (seqs) |s| {
        print("ID:  {s}\nSeq: {s}\n\n", .{s.seqID, s.seq});
    }

    const cname = "copy.fas";
    try bio.writeFastaFile(io, seqs, cname);

    bio.freeFastaArray(gpa, seqs);

    const seqr = try bio.randomFastaArray(gpa, 5);
    print("Sequences in random array\n", .{});
    for (seqr) |s| {
        print("ID:  {s}\nSeq: {s}\n\n", .{s.seqID, s.seq});
    }

    bio.freeFastaArray(gpa, seqr);
}

test "simple test" {
    const gpa = std.testing.allocator;
    var list: std.ArrayList(i32) = .empty;
    defer list.deinit(gpa); // Try commenting this out and see if zig detects the memory leak!
    try list.append(gpa, 42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}

