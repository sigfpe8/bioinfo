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
