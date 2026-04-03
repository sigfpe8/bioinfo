//! This module contains functions that handle protein "motifs".
//! These are a form of (very simplified) regular expression used
//! to describe certain points in amino acid sequences.

const std = @import("std");
const Allocator = std.mem.Allocator;

// The textual language used to describe a motif consists of a string
// of patterns possibly separated by dashes (`-`). The dashes are used
// for human clarity only and are ignored when interpreting the motif.
// Each pattern must match a single amino acid in the test protein.
// 
// The basic pattern is a single uppercase letter that represents one
// particular amino acid. For example:
//      A, E, M or T
// 
// The following patterns can be used in a motif:
// 
// 1. Single uppercase letters describe one particular amino acid
//    (see table below)
// 
//       A   matches only Alanine
//       C   matches only Cysteine
// 
// 2. [Letters...] A sequence of letters inside square brackets describe
//    alternatives. That is, any of those amino acids can match the
//    corresponding amino acid in the test protein.
// 
//       [LVFY] matches leucine, valine, phenylalanine, or tyrosine
// 
// 3. {Letters...} A sequence of letters inside braces describe
//    exclusions. That is, the amino acid in the test protein cannot
//    be any of these.
// 
//        {AC} matches any amino acid except Alanine or Cysteine
// 
// 4. x (lowercase 'X') matches any amino acid
// 
// 5. x(n) matches 'n' amino acids (exact gap)
// 
// 6. x(n,m) matches a minimum of 'n' and a maximum of 'm' amino acids
//
// Amino acids
// Name             Abbr   Code
// ---------------  ----   ----
// Alanine           Ala     A
// Arginine          Arg     R
// Asparagine        Asn     N
// Aspartic acid     Asp     D
// Cysteine          Cys     C
// Glutamic acid     Glu     E
// Glutamine         Gln     Q
// Glycine           Gly     G
// Histidine         His     H
// Isoleucine        Ile     I
// Leucine           Leu     L
// Lysine            Lys     K
// Methionine        Met     M
// Phenylalanine     Phe     F
// Proline           Pro     P
// Serine            Ser     S
// Threonine         Thr     T
// Tryptophan        Trp     W
// Tyrosine          Tyr     Y
// Valine            Val     V
//
// The textual motif is encoded as a binary sequence of bytes in order
// to make the pattern matching algorithm slightly more convenient.
// 
// Basically, the single letter codes for amino acids continue unchanged
// and the alternative and exception patterns are encoded as:
//
//     [ABCD...] -> [nABCD...
//     {ABCD...} -> {nABCD...
//                    \__n__/
//
// where 'n' is the binary number of letters that follow.
//
// The repetition is encoded as
//
//    xnm
//
// where 'n' and 'm' are binary (0-255) versions of the text in `x(n,m)`.

const MotifError = error{
    ExpectedClosingBracket,
    ExpectedClosingBrace,
    EmptyPattern,
};

pub const Motif = struct {
    const Self = @This();

    allocator: Allocator,
    text: []u8,
    pattern: []u8,

    pub fn init(allocator: Allocator) Motif {
        return .{
            .allocator = allocator,
            .text = &.{},
            .pattern = &.{},
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.text);
        self.allocator.free(self.pattern);
    }

    pub fn encode(self: *Self, text: []const u8) !void {
        var buffer: [1024]u8 = undefined;
        self.deinit();  // The motif can be reused with a different pattern
        self.text = try self.allocator.dupe(u8, text);

        var end: usize = 0;
        var inp: usize = 0;
        while (true) {
            if (inp == text.len) {
                break;
            }
            var c = text[inp];
            inp += 1;

            if (c >= 'A' and c <= 'Z') {
                // Single amino acid
                buffer[end] = c;    // Add another single-letter pattern
                end += 1;
                continue;   // Advance to next pattern
            }
            
            if (c == '[') {
                // Alternatives
                buffer[end] = c;
                end += 1;
                const beg = end;    // Remember where to put 'n'
                end += 1;
                while (true) {
                    if (inp == text.len) {
                        return MotifError.ExpectedClosingBracket;
                    }
                    c = text[inp];
                    inp += 1;
                    if (c == ']') {
                        const len: u8 = @truncate(end - beg - 1);
                        if (len == 0) {
                            return MotifError.EmptyPattern;
                        }
                        buffer[beg] = len;
                        break;
                    }
                    buffer[end] = c;    // Add another letter to the alternative set
                    end += 1;
                }
                continue;   // Advance to next pattern
            }

            if (c == '{') {
                // Exceptions
                buffer[end] = c;
                end += 1;
                const beg = end;    // Remember where to put 'n'
                end += 1;
                while (true) {
                    if (inp == text.len) {
                        return MotifError.ExpectedClosingBrace;
                    }
                    c = text[inp];
                    inp += 1;
                    if (c == '}') {
                        const len: u8 = @truncate(end - beg - 1);
                        if (len == 0) {
                            return MotifError.EmptyPattern;
                        }
                        buffer[beg] = len;
                        break;
                    }
                    buffer[end] = c;    // Add another letter to the exclusion set
                    end += 1;
                }
                continue;   // Advance to next pattern
            }
        }

        self.pattern = try self.allocator.dupe(u8, buffer[0..end]);
    }

    pub fn decode(self: *Self) void {
        var i: usize = 0;

        while (i < self.pattern.len) {
            var c = self.pattern[i];
            i += 1;
            if (c >= 'A' and c <= 'Z') {
                std.debug.print("{c}", .{c});
                continue;
            }

            // We assume the pattern has been correctly encoded so
            // that we don't need to check whether i > pattern.len 
            if (c == '[') {
                var n = self.pattern[i];
                i += 1;
                std.debug.print("[{d}:", .{n});
                while (n > 0) : (n -= 1) {
                    c = self.pattern[i];
                    i += 1;
                    std.debug.print("{c}", .{c});
                }
                std.debug.print("]", .{});
                continue;
            }

            if (c == '{') {
                var n = self.pattern[i];
                i += 1;
                std.debug.print("{{{d}:", .{n});
                while (n > 0) : (n -= 1) {
                    c = self.pattern[i];
                    i += 1;
                    std.debug.print("{c}", .{c});
                }
                std.debug.print("}}", .{});
                continue;
            }
        }
        std.debug.print("\n", .{});
    }

    pub fn match(self: *Self, text: []const u8) bool {
        var it: usize = 0;  // Text to match
        var ip: usize = 0;  // Pattern

        while (true) {
            if (ip == self.pattern.len) {
                return true;    // Exhausted the pattern    
            }
            const p = self.pattern[ip];
            ip += 1;
            if (it == text.len) {
                return false;   // Exhausted the text but not the pattern
            }
            const t = text[it];
            it += 1;
            if (p >= 'A' and p <= 'Z') {
                if (p == t)
                    continue;   // Match, advance to the next pattern
                return false;
            }
            if (p == '[') {
                // [ n A B ...
                //     \_ n _/
                const n = self.pattern[ip];
                ip += 1;
                var i: usize = ip;
                ip += n;    // Next pattern
                while (i < ip) : (i += 1) {
                    if (self.pattern[i] == t) { // Found an alternative matching
                        break;
                    }
                } else {    // Exhausted the alternatives without any match
                    return false;
                }
                continue;   // Match, advance to the next pattern
            }
            if (p == '{') {
                // { n A B ...
                //     \_ n _/
                const n = self.pattern[ip];
                ip += 1;
                var i: usize = ip;
                ip += n;    // Next pattern
                while (i < ip) : (i += 1) {
                    if (self.pattern[i] == t) { // Matched a forbidden letter
                        return false;
                    }
                }
                continue;   // Match, advance to the next pattern
            }
        }
    }
};

