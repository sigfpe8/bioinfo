const std = @import("std");
const bio = @import("bioinfo");

const Io = std.Io;
const Reader = Io.Reader;
const Writer = Io.Writer;
const Allocator = std.mem.Allocator;
const print = std.debug.print;

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;
    const seq = try bio.randomSeq(gpa, 50);
    defer gpa.free(seq);
    std.debug.print("Random sequence: {s}\n", .{seq});

    const fname = "datasets/example.fas";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    bio.printLines(lines);
    // const seqs = try bio.readFastaFile(io, gpa, fname);

    // print("Sequences in file '{s}':\n", .{fname});
    // for (seqs) |s| {
    //     print("ID:  {s}\nSeq: {s}\n\n", .{s.seqID, s.seq});
    // }

    // const cname = "copy.fas";
    // try bio.writeFastaFile(io, seqs, cname);

    // bio.freeFastaArray(gpa, seqs);

    // const seqr = try bio.randomFastaArray(gpa, 5);
    // print("Sequences in random array\n", .{});
    // for (seqr) |s| {
    //     print("ID:  {s}\nSeq: {s}\n\n", .{s.seqID, s.seq});
    // }

    // bio.freeFastaArray(gpa, seqr);
}

