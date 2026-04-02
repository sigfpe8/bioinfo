const std = @import("std");
const Io = std.Io;
const Reader = Io.Reader;
const Writer = Io.Writer;
const Allocator = std.mem.Allocator;
const print = std.debug.print;
const swap = std.mem.swap;
const assert = std.debug.assert;

/// Utilities
/// zig 0.16.0

/// Binomial coefficient
/// n choose k
pub fn binomial(n: usize, k: usize) usize {
    // (n * (n - 1) * ... * (n - k + 1)) / (k * (k - 1) * ... * 1)
    if (n == k or k == 0) {
        return 1;
    }

    assert(n > k);

    var t: usize = n - 1;
    var num: usize = n;

    // Numerator = n * (n - 1) * ... * (n - k + 1)
    // k factors, (k-1) multiplications
    while (t > (n - k)) : (t -= 1) {
        num *= t;
    }

    // Denominator =  k * (k - 1) * ... * 1
    // k factors, (k-1) multiplications
    var den: usize = k;
    t = k - 1;
    while (t > 1) : (t -= 1) {
        den *= t;
    }

    return num / den;
}

/// Print Pascal's Triangle of order n
pub fn pascalsTriang(n: usize) void {
    for (0..n+1) |i| {
        for (0..i+1) |j| {
            print("{d:3} ", .{ binomial(i,j)});
        }
        print("\n", .{});
    }
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
pub fn printLines(lines: [][]u8) void {
    for (lines) |line| {
        print("{s}\n", .{line});
    }
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
            // const file = try Io.Dir.cwd().openFile(io, fname, .{ .mode = .read_only });
            // defer file.close(io);

            // // Buffer to hold 1 line
            // const buffer: []u8 = try self.allocator.alloc(u8, 1024);
            // defer self.allocator.free(buffer);

            // var file_reader: Io.File.Reader = .init(file, io, buffer);
            // var reader: *Reader = &file_reader.interface;

            // const line = try reader.takeDelimiterInclusive('\n');

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
        pub fn printInts(self: *Self) void {
            for (self.items) |int| {
                print("{d}\n", .{int});
            }
        }
    };
}

/// Permutation iterator
/// Heap's algorithm as implemented by @DRcrazy3
/// https://youtu.be/atD5EWada7k?si=nci2lVksjs2AAcAc
pub fn PermIterator(comptime T: type) type {
    return struct {
        const Self = @This();

        items: []T,
        counter: []usize, // c
        ptr: usize, // i
        first: bool,

        pub fn init(alloc: Allocator, slice: []const T) !Self {
            const items = try alloc.dupe(T, slice);
            const counter = try alloc.alloc(usize, slice.len);
            @memset(counter, 0);

            return .{
                .items = items,
                .counter = counter,
                .ptr = 1,
                .first = true,
            };
        }

        pub fn deinit(self: *Self, alloc: Allocator) void {
            alloc.free(self.items);
            alloc.free(self.counter);
        }

        pub fn next(self: *Self) ?[]T {
            if (self.first) {
                self.first = false;
                return self.items;
            }

            while (self.ptr < self.items.len) {
                const count = self.counter[self.ptr];

                if (count < self.ptr) {
                    // swap or rotate
                    const j = if (self.ptr % 2 == 0) 0 else count;
                    swap(T, &self.items[self.ptr], &self.items[j]);

                    self.counter[self.ptr] += 1;
                    self.ptr = 1;
                    return self.items;
                }

                self.counter[self.ptr] = 0;
                self.ptr += 1;
            }

            return null;
        }
    };
}
