//! Support for FASTA files

const std = @import("std");

const Io = std.Io;
const Reader = Io.Reader;
const Writer = Io.Writer;
const Allocator = std.mem.Allocator;

const FastaError = error{
        FailedToLoadProtein,
};

// This struct represents a single sequence in a FASTA file.
// A file can contain multiple sequences.
pub const Fasta = struct {
    seqID: []u8,
    seq: []u8,
};

// pub fn randomFastaArray(allocator: Allocator, n: usize) ![]Fasta {
//     var prng = std.Random.DefaultPrng.init(0); 
//     const rand = prng.random();

//     var list: std.ArrayList(Fasta) = .empty;
//     var i: usize = 1;

//     while (i <= n) : (i += 1) {
//         const seqId = try std.fmt.allocPrint(allocator, "Sequence {d}", .{i});
//         const len = rand.uintLessThan(usize, 100) + 50;
//         const seq = try randomSeq(allocator, len);
//         try list.append(allocator, .{ .seqID = seqId, .seq = seq });
//     }

//     return list.toOwnedSlice(allocator);
// }

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
                return FastaError.FailedToLoadProtein;
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
fn fetchProtein(io: Io, gpa: Allocator, prot_name: []const u8) ![]u8 {
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
        return FastaError.FailedToLoadProtein;
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

/// Read the first line from a file
pub fn readLine(io: Io, allocator: Allocator, fname: []const u8) ![]u8 {
    var buffer: [1024]u8 = undefined;
    const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
    defer file.close(io);

    var file_reader: Io.File.Reader = .init(file, io, &buffer);
    var reader: *Reader = &file_reader.interface;
    const line = try reader.takeDelimiterInclusive('\n');

    return allocator.dupe(u8, line[0..line.len - 1]);   // Trim \n
}

/// Read all lines in a file into a slice of text slices
pub fn readLines(io: Io, allocator: Allocator, fname: []const u8) ![][]u8 {
    const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
    defer file.close(io);

    // Buffer to hold 1 line, which can be very large
    const buffer: []u8 = try allocator.alloc(u8, 16*1024);
    defer allocator.free(buffer);

    var file_reader: Io.File.Reader = .init(file, io, buffer);
    var reader: *Reader = &file_reader.interface;

    var lines: std.ArrayList([]u8) = .empty;
    var line: []u8 = undefined;

    while (true) {
        line = reader.takeDelimiterInclusive('\n') catch |err| {
            switch (err) {
                error.EndOfStream => break,
                else => return err,
            }
        };
        // Add one more line
        line = try allocator.dupe(u8, line[0..line.len - 1]); // Trim \n
        try lines.append(allocator, line);
    }

    return lines.toOwnedSlice(allocator);
}

/// Free lines
pub fn freeLines(allocator: Allocator, lines: [][]u8) void {
    for (lines) |line| {
        allocator.free(line);
    }
    allocator.free(lines);
}

/// Print lines
pub fn printLines(lines: [][]u8, writer: *Writer) !void {
    for (lines) |line| {
        try writer.print("{s}\n", .{line});
    }
    try writer.flush();
}

pub fn IntsReader(comptime T: type) type {
    return struct {
        const Self = @This();

        allocator: Allocator,
        items: []T,

        pub fn init(allocator: Allocator) Self {
            return .{
                .allocator = allocator,
                .items = &.{},
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.items);
            self.items = &.{};
        }

        /// Read one line from a file and return a list of the parsed integers in it
        pub fn readFile(self: *Self, io: Io, fname: []const u8) !void {
            const line = try readLine(io, self.allocator, fname);
            defer self.allocator.free(line);

            try self.readBuffer(line);
        }

        pub fn readBuffer(self: *Self, buffer: []const u8) !void {
            var ints: std.ArrayList(T) = .empty;
            var beg: usize = 0;
            var end: usize = 0;

            // Parse buffer, e.g.: " 27   35   1025  "
            while (beg < buffer.len) : (beg = end) {
                // Skip whitespace at the beginning
                while (beg < buffer.len) {
                    if (std.ascii.isWhitespace(buffer[beg])) {
                        beg += 1;
                    } else {
                        break;
                    }
                }
                // End of the buffer?
                if (beg == buffer.len) {
                    break;  // Yes, no more ints to parse
                }

                // Find end of number
                end = beg + 1;  // Ok if end == buffer.len
                while (end < buffer.len) {
                    if (std.ascii.isDigit(buffer[end])) {
                        end += 1;
                    } else {
                        break;
                    }
                }
                const int = try std.fmt.parseInt(T, buffer[beg..end], 10);
                try ints.append(self.allocator, int);
            }

            self.items = try ints.toOwnedSlice(self.allocator);
        }

        /// Print integers
        pub fn printInts(self: *Self, writer: *Writer) !void {
            for (self.items) |int| {
                try writer.print("{d}\n", .{int});
            }
            try writer.flush();
        }
    };
}