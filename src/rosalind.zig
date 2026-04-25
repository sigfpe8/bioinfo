//! Rosalind problems and solutions
//! Bioinformatics stronghold: https://rosalind.info/problems/list-view/

const std = @import("std");
const bio = @import("root.zig");
const Fasta = bio.Fasta;
const binomial = bio.binomial;
const aacTable = bio.aacTable;
const stopCodons = bio.stopCodons;

const Io = std.Io;
const File = Io.File;
const Reader = Io.Reader;
const Writer = Io.Writer;
const StringSet = bio.StringSet;

const Allocator = std.mem.Allocator;
const assert = std.debug.assert;

var gpa: Allocator = undefined;
var io: std.Io = undefined;
var stdout_buffer: [1024]u8 = undefined;
var file_writer: File.Writer = undefined;
var stdout: *Writer = undefined; 

pub fn main(init: std.process.Init) !void {
    gpa = init.gpa;
    io = init.io;
    file_writer = std.Io.File.Writer.init(.stdout(), io, &stdout_buffer);
    stdout = &file_writer.interface;

    // try problemDNA();
    // try problemRNA();
    // try problemREVC();
    // try problemGC();
    // try problemFIB();
    // try problemHAMM();
    // try problemPERM();
    // try problemIPRB();
    // try problemSUBS();
    // try problemPROT();
    // try problemPROT2();
    // try problemREVP();
    // try problemCONS();
    // try problemMPRT();
    // try problemGRPH();
    // try problemMRNA();
    // try problemFIBD();
    // try problemPRTM();
    // try problemSPLC();
    // try problemORF();
    // try problemLEXF();
    // try problemLCSM();
    // try problemIEV();
    try problemLIA();
}

/// Simple print() to stdout ignoring errors
fn print(comptime fmt: []const u8, args: anytype) void {
    stdout.print(fmt, args) catch {};
    stdout.flush() catch {};
}

// Each problemXYZ() function below gathers the necessary data for
// problem XYZ and then calls the corresponding solveXYZ() function
// to actually solve the given problem. All these functions assume
// that the datasets contain valid data and are not empty. If a more
// robust program is desired, one can test these and other conditions
// before calling solveXYZ(). For example, after reading a FASTA file:
//
//    if (seqs.len > 0) {
//        solveXYZ(seqs);   // Expects a non-empty slice
//    }

fn problemDNA() !void {
    const fname = "datasets/rosalind_dna.txt";
    const line = try bio.readLine(io, gpa, fname);
    defer gpa.free(line);

    solveDNA(line);
}

fn problemRNA() !void {
    const fname = "datasets/rosalind_rna.txt";
    const line = try bio.readLine(io, gpa, fname);
    defer gpa.free(line);

    try solveRNA(line);
}

fn problemREVC() !void {
    const fname = "datasets/rosalind_revc.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    try solveREVC(lines[0]);
}

fn problemGC() !void {
    const fname = "datasets/rosalind_gc.txt";
    const seqs = try bio.readFastaFile(io, gpa, fname);
    defer bio.freeFastaArray(gpa, seqs);

    try solveGC(seqs);
}

fn problemFIB() !void {
    const fname = "datasets/rosalind_fib.txt";
    var ints = bio.IntsReader(usize).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    solveFIB(ints.items[0], ints.items[1]);
}

fn problemHAMM() !void {
    const fname = "datasets/rosalind_hamm.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    solveHAMM(lines[0], lines[1]);
}

fn problemPERM() !void {
    const fname = "datasets/rosalind_perm.txt";
    var ints = bio.IntsReader(usize).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    try solvePERM(ints.items[0]);
}

fn problemIPRB() !void {
    const fname = "datasets/rosalind_iprb.txt";
    var ints = bio.IntsReader(usize).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    solveIPRB(ints.items[0], ints.items[1], ints.items[2]);
}

fn problemSUBS() !void {
    const fname = "datasets/rosalind_subs.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    solveSUBS(lines[0], lines[1]);
}

fn problemPROT() !void {
    const fname = "datasets/rosalind_prot.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    try solvePROT(lines[0]);
}

fn problemPROT2() !void {
    const fname = "datasets/sample_prot.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    try solvePROT2(lines[0]);
}

fn problemREVP() !void {
    const fname = "datasets/rosalind_revp.txt";
    const seqs = try bio.readFastaFile(io, gpa, fname);
    defer bio.freeFastaArray(gpa, seqs);

    solveREVP(seqs[0].seq);
}

fn problemCONS() !void {
    const fname = "datasets/rosalind_cons.txt";
    const seqs = try bio.readFastaFile(io, gpa, fname);
    defer bio.freeFastaArray(gpa, seqs);

    const seqs2 = try bio.fastaToSequences(gpa, seqs);
    defer bio.freeLines(gpa, seqs2);

    try solveCONS(seqs2);
}

fn problemMPRT() !void {
    const fname = "datasets/rosalind_mprt.txt";
    const list = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, list);

    try solveMPRT(list);
}

fn problemGRPH() !void {
    const fname = "datasets/rosalind_grph.txt";
    const seqs = try bio.readFastaFile(io, gpa, fname);
    defer bio.freeFastaArray(gpa, seqs);

    solveGRPH(seqs, 3);
}

fn problemMRNA() !void {
    const fname = "datasets/rosalind_mrna.txt";
    const prot = try bio.readLine(io, gpa, fname);
    defer gpa.free(prot);

    solveMRNA(prot);
}

fn problemFIBD() !void {
    const fname = "datasets/rosalind_fibd.txt";
    var ints = bio.IntsReader(usize).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    solveFIBD(ints.items[0], ints.items[1]);
}

fn problemPRTM() !void {
    const fname = "datasets/rosalind_prtm.txt";
    const prot = try bio.readLine(io, gpa, fname);
    defer gpa.free(prot);

    solvePRTM(prot);
}

fn problemSPLC() !void {
    const fname = "datasets/rosalind_splc.txt";
    const seqs = try bio.readFastaNoIdFile(io, gpa, fname);
    defer bio.freeLines(gpa, seqs);

    try solveSPLC(seqs[0], seqs[1..]);
}

fn problemORF() !void {
    const fname = "datasets/rosalind_orf.txt";
    const seqs = try bio.readFastaNoIdFile(io, gpa, fname);
    defer bio.freeLines(gpa, seqs);

    try solveORF(seqs[0]);
}

fn problemLEXF() !void {
    const fname = "datasets/rosalind_lexf.txt";
    const lines = try bio.readLines(io, gpa, fname);
    defer bio.freeLines(gpa, lines);

    // Assume there's a space between the symbols, eg: "A B C D"
    const len = (lines[0].len + 1) / 2;

    var alphabet = try gpa.alloc(u8, len);
    defer gpa.free(alphabet);

    // Create alphabet without the spaces
    var i: usize = 0;
    for (lines[0]) |c| {
        if (c != ' ') {
            alphabet[i] = c;
            i += 1;
        }
    }

    const n: usize = try std.fmt.parseInt(usize, lines[1], 10);

    try solveLEXF(alphabet, n);
}

fn problemLCSM() !void {
    const fname = "datasets/rosalind_lcsm.txt";
    const seqs = try bio.readFastaNoIdFile(io, gpa, fname);
    defer bio.freeLines(gpa, seqs);

    try solveLCSM(seqs);
}

fn problemIEV() !void {
    const fname = "datasets/rosalind_iev.txt";
    var ints = bio.IntsReader(i32).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    try solveIEV(ints.items);
}

fn problemLIA() !void {
    const fname = "datasets/rosalind_lia.txt";
    var ints = bio.IntsReader(usize).init(gpa);
    defer ints.deinit();
    try ints.readFile(io, fname);

    solveLIA(ints.items[0], ints.items[1]);
}
// ---------------------------------------------------------

/// DNA - Counting DNA Nucleotides
pub fn solveDNA(seq: []const u8) void {
    const nA, const nC, const nG, const nT = bio.countBases(seq);
    print("{d} {d} {d} {d}\n", .{nA,nC,nG,nT});
}

/// RNA - Transcribing DNA into RNA
pub fn solveRNA(seq: []const u8) !void {
    const rna = try bio.transcribe(gpa, seq);
    print("{s}\n", .{rna});
    gpa.free(rna);
}

/// REVC - Complementing a Strand of DNA
fn solveREVC(seq: []const u8) !void {
    const revC = try bio.revComplement(gpa, seq);
    print("{s}\n", .{revC});
    gpa.free(revC);
}

/// GC - Computing GC content
fn solveGC(seqs: []Fasta) !void {
    var gc = bio.gcContent(seqs[0].seq);
    var id: usize = 0;

    for (seqs[1..], 1..) |s, i| {
        const x = bio.gcContent(s.seq);
        if (x > gc) {
            gc = x;
            id = i;
        }
    }
    print("{s}\n", .{seqs[id].seqID});
    print("{d:.6}\n", .{gc});
}

/// FIB - Rabbits andRecurrence Relations
pub fn solveFIB(n: usize, k: usize) void {
    var fn2: usize = 1;
    var fn1: usize = 1;
    var f: usize = 1;
    var i: usize = 3;

    while (i <= n) : (i += 1) {
        // Recurrence relation: Fn = Fn-1 + k * Fn-2
        f = fn1 + k * fn2;
        fn2 = fn1;
        fn1 = f;
    }

    print("{d}\n", .{f});
}

/// HAMM - Counting Point Mutations
fn solveHAMM(seq1: []const u8, seq2: []const u8) void {
    const hamd = bio.hammingDist(seq1, seq2);

    print("{d}\n", .{hamd});
}

/// PERM - Enumerating Gene Orders
pub fn solvePERM(n: usize) !void {
    var perm1 = try gpa.alloc(u8, n);
    defer gpa.free(perm1);

    var fact: usize = 1;
    for (0..n) |i| {
        const j = i + 1;
        perm1[i] = @truncate(j);
        fact = fact * j;
    }
    print("{d}\n", .{fact});

    var pit = try bio.PermIterator(u8).init(gpa, perm1);
    defer pit.deinit(gpa);

    while (pit.next()) |perm| {
        for (perm) |p| {
            print("{d} ", .{p});
        }
        print("\n", .{});
    }
}   

/// IPRB - Mendel's First Law
pub fn solveIPRB(k: usize, m: usize, n:usize) void {
    // Number of possible pairs = (k + m + n) choose 2
    const t = binomial(k+m+n,2);
    const tf: f64 = @floatFromInt(t);

    // Number of pairs with at least 1 individual from k (homozygous dominant)
    const a = binomial(k,2) + k * m + k * n;
    const af: f64 = @floatFromInt(a);

    // Number of pairs with two individuals from m (heterozygous)
    const b = binomial(m,2);
    const bf: f64 = @floatFromInt(b);

    // Number of pairs with 1 individual from m and the other from n
    const c = m * n;
    const cf: f64 = @floatFromInt(c);

    const p = (af/tf) + (bf/tf)*0.75 + (cf/tf)*0.5;
    print("{d:.5}\n", .{p});
}

/// SUBS - Finding a Motif in DNA 
pub fn solveSUBS(seq: []const u8, sub: []const u8) void {
    var base: usize = 0;
    while (true) {
        const pos = std.mem.indexOf(u8, seq[base..], sub);
        if (pos) |p| {
            print("{d} ", .{base+p+1});
            base += p + 1;  // Restart at the next character
        } else {
            print("\n", .{});
            break;
        }
    }
}

/// PROT - Translating RNA into Protein
fn solvePROT(rna: []const u8) !void {
    assert(rna.len % 3 == 0);

    var map = try bio.makeRNACodonMap(gpa);
    defer map.deinit();

    var prot = try gpa.alloc(u8, rna.len / 3);
    defer gpa.free(prot);

    var i: usize = 0;
    while (i < rna.len) : (i += 3) {
        const cod = rna[i..i+3];        // The 3-letter codon
        const aac = map.get(cod).?;     // Corresponding amino acid
        prot[i / 3] = aac;
    }

    if (prot[prot.len-1] == '.') {
        print("{s}\n", .{prot[0..prot.len-1]}); 
    } else {
        print("{s}\n", .{prot}); 
    }
}

fn solvePROT2(rna: []const u8) !void {
    const prot = try bio.translateRNA(gpa, rna);
    defer gpa.free(prot);
    
    print("{s}\n", .{prot});
}

/// REVP - Locating Restriction Sites
fn solveREVP(seq: []const u8) void {
    var base: usize = 0;

    while (base < seq.len) : (base += 1) {
        var i: usize = 4;
        while (i <= 12) : (i += 2) {
            if (base + i > seq.len) {
                break;
            }
            if (bio.isRevPalindrome(seq[base..base+i])) {
                print("{d} {d}\n", .{base+1, i});
            }
        }
    }
}

/// CONS - Consensus and Profile
fn solveCONS(seqs: [][]u8) !void {
    const prof = try bio.profile(gpa, seqs);
    defer bio.freeProfile(gpa, prof);

    const cons = try bio.consensus(gpa, prof);
    defer gpa.free(cons);

    // Print the consensus string
    print("{s}\n", .{cons});

    // Print the profile matrix
    for (prof, 0..) |cnts, i| {
        const sym: u8 = switch (i) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            else => unreachable,
        };

        print("{c}:", .{sym});
        for (cnts) |cnt| {
            print(" {d}", .{cnt});
        }
        print("\n", .{});
    }
}

/// MPRT - Finding a Protein Motif
fn solveMPRT(list: [][] u8) !void {
    var motif = bio.Motif.init(gpa);
    defer motif.deinit();
    // N-glycosylation motif
    try motif.encode("N{P}[ST]{P}");

    for (list) |prot_name| {
        const fasta = try bio.getProtein(io, gpa, prot_name);
        const seq = fasta[0].seq;
        const locs  = try motif.locations(seq);

        if (locs.len > 0) {
            print("{s}\n", .{prot_name});
            for (locs) |loc| {
                print("{d} ", .{loc});
            }
            print("\n", .{});
        }
        // Free resources for this protein
        gpa.free(locs);
        bio.freeFastaArray(gpa, fasta);
    }
}

/// GRPH - Overlap Graphs
fn solveGRPH(seqs: []Fasta, k: usize) void {
    var i: usize = 0;
    var j: usize = 0;

    while (i < seqs.len - 1) : (i += 1) {
        j = i + 1;
        const s1 = seqs[i].seq;
        while (j < seqs.len) : (j += 1) {
            const s2 = seqs[j].seq;
            // s1 -> s2 ?
            if (std.mem.eql(u8, s1[s1.len - k..], s2[0..k])) {
                print("{s} {s}\n", .{seqs[i].seqID, seqs[j].seqID});
            }
            // s2 -> s1 ?
            if (std.mem.eql(u8, s2[s2.len - k..], s1[0..k])) {
                print("{s} {s}\n", .{seqs[j].seqID, seqs[i].seqID});
            }
        }
    }
}

/// MRNA - Inferring mRNA from Protein
fn solveMRNA(prot: []const u8) void {
    // Assume the protein string contains at least one amino acid
    // Assume the string contains only letters A-Z, no stop "." at the end
    assert(prot.len > 0);
    var count: usize = stopCodons;   // 3

    for (prot) |a| {
        count = (count * aacTable[a - 'A']) % 1_000_000;
    }

    print("{d}\n", .{count});
}

/// FIBD - Mortal Fibonacci Rabbits
fn solveFIBD(n: usize, m: usize) void {
    assert(m <= 20);
    print("n={d}, m={d}\n", .{n,m});
    var age = [_]usize{ 0 } ** 20;
    age[0] = 1;
    var i: usize = 1;
    while (i < n) : (i += 1) {
        // Each generation at least 1 month old
        // generates new children
        const nb = allGens(age[1..m]);
        // Each generation gets 1 month older
        // Oldest generation dies (m-month old)
        var j: usize = m - 1;
        while (j > 0) : (j -= 1) {
            age[j] = age[j - 1];
        }
        // Newborns
        age[0] = nb;
    }
    const count = allGens(age[0..m]);
    print("{d}\n", .{count});
 }

// Sums number of rabitts for all generations
fn allGens(age: []usize) usize {
    var count: usize = 0;
    for (age) |c| {
        count += c;
    }
    return count;
}

/// PRTM - Calculating Protein Mass
fn solvePRTM(pro: []const u8) void {
    const mass = bio.proteinMass(pro);
    print("{d:.3}\n", .{mass});
}

/// SPLC - RNA Splicing
fn solveSPLC(dna: []const u8, introns: [][]u8) !void {
    const prot = try bio.spliceDNA(gpa, dna, introns);
    defer gpa.free(prot);

    print("{s}\n", .{prot});
}

/// ORF - Open Reading Frames
fn solveORF(dna: []const u8) !void {
    var set = try bio.findORFs(gpa, dna);
    defer set.deinit();

    var it = set.iterator();
    while (it.next()) |key_ptr| {
        print("{s}\n", .{key_ptr.*});
    }
}

/// LEXF - Enumerating k-mers Lexicographically
fn solveLEXF(alphabet: []const u8, n: usize) !void {
    const string = try gpa.alloc(u8, n);
    defer gpa.free(string);

    var sym: u8 = 0;    // Symbol index from alphabet [0..alphabet.len]
    var ind: usize = 0; // Index in string [0..n]

    next_string:
    while (true) {
        // Fill following positions with first symbol (0)
        while (ind < n) : (ind += 1) {
            string[ind] = 0;
        }

        // Print the string
        for (string) |i| {
            print("{c}", .{alphabet[i]});
        }
        print("\n", .{});

        // Bump last symbol, backtracking as necessary
        // Enter with ind = n
        while (true) {
            // Next symbol
            ind -= 1;
            sym = string[ind] + 1;
            if (sym == alphabet.len) {
                // Already highest symbol, must backtrack
                if (ind == 0) {
                    // Reached first column, cannot backtrack
                    break :next_string;
                }
            } else {
                break;
            }
        }

        string[ind] = sym;
        ind += 1;
    }
}

/// Return true if all the strings contain the substring
fn allContain(strs: [][]u8, subs: []const u8) bool {
    for (strs) |s| {
        if (std.mem.find(u8, s, subs) == null) {
            return false;
        }
    }

    return true;
}

/// LCSM - Finding a Shared Motif
fn solveLCSM(seqs: [][]u8) !void {
    // Find the shortest sequence
    // The longest common substring cannot be longer than it
    var len = seqs[0].len;
    var ind: usize = 0;

    for (seqs[1..], 1..) |s, i| {
        if (s.len < len) {
            len = s.len;
            ind = i;
        }
    }

    const needle = seqs[ind];

    // print("len={d}, ind={d}, needle={s}\n", .{len,ind,needle});

    while (len > 0) : (len -= 1) {
        var base: usize = 0;
        while (base + len <= needle.len) : (base += 1) {
            if (allContain(seqs, needle[base..base+len])) {
                // Found one of the longest common substrings
                print("{s}\n", .{needle[base..base+len]});
                return;
            }
        }

    } else {
        print("The strings have no common substring!\n", .{});
        return;
    }
}

/// IEV - Calculating Expected Offspring
fn solveIEV(ints: []i32) !void {
    assert(ints.len == 6);

    const prob = [_]f64{1.0, 1.0, 1.0, 0.75, 0.5, 0};

    var exp: f64 = 0.0;

    for (ints, 0..) |v, i| {
        const fv: f64 = @floatFromInt(v);
        // print("{d} - {d}\n", .{fv, prob[i]});
        exp += 2.0 * fv * prob[i];
    }

    print("{d}\n", .{exp});
}

/// LIA - Independent Alleles
fn solveLIA(k: usize, n: usize) void {
    var prob: f64 = 0.0;

    // Number of offsprings after k generations (2**k)
    const offs: usize = std.math.pow(usize, 2, k);

    var i: usize = n;
    while (i <= offs) : (i += 1) {
        const t1: f64 = @floatFromInt(bio.nChooseK(offs, i));
        const t2: f64 = std.math.pow(f64, 0.25, @as(f64, @floatFromInt(i)));
        const t3: f64 = std.math.pow(f64, 0.75, @as(f64, @floatFromInt(offs - i)));
        prob += t1 * t2 * t3;
    }

    print("{d:.3}\n", .{prob});
}