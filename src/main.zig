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

    const fname = "datasets/sample_grph.txt";
    const seqs = try bio.readFastaNoIdFile(io, gpa, fname);
    defer bio.freeLines(gpa, seqs);

    print("Sequences in file '{s}':\n", .{fname});
    for (seqs) |s| {
        print("Seq: {s}\n", .{s});
    }

    // var motif = bio.Motif.init(gpa);
    // defer motif.deinit();
    // try motif.encode("ABC[DEF]G{HIJKL}");
    // // try motif.encode("A{GFDK}C");
    // const decod = try motif.decode();
    // defer gpa.free(decod);
    // print("Motif = {s}\n", .{decod});
    // const text = "ABCEGM";
    // print("Text  = {s}\n", .{text});
    // if (motif.match(text)) {
    //     print("Success\n", .{});
    // } else {
    //     print("Failure\n", .{});        
    // }
}

