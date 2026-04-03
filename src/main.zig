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

    const fname = "datasets/example.fas";
    const seqs = try bio.readFastaFile(io, gpa, fname);
    defer bio.freeFastaArray(gpa, seqs);

    print("Sequences in file '{s}':\n", .{fname});
    for (seqs) |s| {
        print("ID:  {s}\nSeq: {s}\n\n", .{s.seqID, s.seq});
    }
}

