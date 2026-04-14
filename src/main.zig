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

    const fname = "datasets/sample_splc.txt";
    const seqs = try bio.readFastaNoIdFile(io, gpa, fname);
    defer bio.freeLines(gpa, seqs);

    for (seqs) |seq| {
        print("{s}\n", .{seq});
    }
}

