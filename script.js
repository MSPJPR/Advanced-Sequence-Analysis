// Global Variables
const gapPenalty = -1; // Gap penalty for alignment

// Event Listeners
document.getElementById('processSequences').addEventListener('click', processSequences);
document.getElementById('alignSequences').addEventListener('click', alignSequences);
document.getElementById('analyzeCodonUsage').addEventListener('click', analyzeCodonUsage);
document.getElementById('searchMotif').addEventListener('click', searchMotif);

// Helper Functions
function readFile(file, callback) {
    const reader = new FileReader();
    reader.onload = () => callback(reader.result);
    reader.readAsText(file);
}

// Upload and Display Sequences
function processSequences() {
    const fileInput = document.getElementById('sequenceFile').files[0];
    if (fileInput) {
        readFile(fileInput, (text) => {
            document.getElementById('sequenceInput').value = text;
        });
    }
}

// Sequence Alignment (Global and Local)
function alignSequences() {
    const seq1 = document.getElementById('alignmentSeq1').value.trim();
    const seq2 = document.getElementById('alignmentSeq2').value.trim();
    const alignmentType = document.getElementById('alignmentType').value;
    let result = '';

    if (seq1 && seq2) {
        if (alignmentType === 'global') {
            result = needlemanWunsch(seq1, seq2);
        } else if (alignmentType === 'local') {
            result = smithWaterman(seq1, seq2);
        }
    } else {
        result = 'Please enter both sequences.';
    }

    document.getElementById('alignmentResult').innerHTML = `<pre>${result}</pre>`;
}

// Needleman-Wunsch Algorithm (Global Alignment)
function needlemanWunsch(seq1, seq2) {
    const n = seq1.length + 1;
    const m = seq2.length + 1;
    const dp = Array.from({ length: n }, () => Array(m).fill(0));
    const traceback = Array.from({ length: n }, () => Array(m).fill(''));

    // Initialize scoring matrix
    for (let i = 0; i < n; i++) dp[i][0] = i * gapPenalty;
    for (let j = 0; j < m; j++) dp[0][j] = j * gapPenalty;

    // Fill scoring matrix
    for (let i = 1; i < n; i++) {
        for (let j = 1; j < m; j++) {
            const match = dp[i - 1][j - 1] + (seq1[i - 1] === seq2[j - 1] ? 1 : -1);
            const deleteGap = dp[i - 1][j] + gapPenalty;
            const insertGap = dp[i][j - 1] + gapPenalty;
            dp[i][j] = Math.max(match, deleteGap, insertGap);

            // Traceback
            if (dp[i][j] === match) traceback[i][j] = '↖';
            else if (dp[i][j] === deleteGap) traceback[i][j] = '↑';
            else traceback[i][j] = '←';
        }
    }

    return formatAlignment(dp, traceback, seq1, seq2, n, m);
}

// Smith-Waterman Algorithm (Local Alignment)
function smithWaterman(seq1, seq2) {
    const n = seq1.length + 1;
    const m = seq2.length + 1;
    const dp = Array.from({ length: n }, () => Array(m).fill(0));
    const traceback = Array.from({ length: n }, () => Array(m).fill(''));

    // Fill scoring matrix
    let maxScore = 0;
    let maxPos = [0, 0];
    for (let i = 1; i < n; i++) {
        for (let j = 1; j < m; j++) {
            const match = dp[i - 1][j - 1] + (seq1[i - 1] === seq2[j - 1] ? 1 : -1);
            const deleteGap = dp[i - 1][j] + gapPenalty;
            const insertGap = dp[i][j - 1] + gapPenalty;
            dp[i][j] = Math.max(0, match, deleteGap, insertGap);

            // Traceback
            if (dp[i][j] === match) traceback[i][j] = '↖';
            else if (dp[i][j] === deleteGap) traceback[i][j] = '↑';
            else if (dp[i][j] === insertGap) traceback[i][j] = '←';

            // Track max score
            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxPos = [i, j];
            }
        }
    }

    return formatAlignment(dp, traceback, seq1, seq2, n, m, maxPos);
}

// Format Alignment for Display
function formatAlignment(dp, traceback, seq1, seq2, n, m, maxPos = [n - 1, m - 1]) {
    let alignment1 = '';
    let alignment2 = '';
    let i = maxPos[0];
    let j = maxPos[1];

    while (i > 0 || j > 0) {
        if (traceback[i][j] === '↖') {
            alignment1 = seq1[i - 1] + alignment1;
            alignment2 = seq2[j - 1] + alignment2;
            i--;
            j--;
        } else if (traceback[i][j] === '↑') {
            alignment1 = seq1[i - 1] + alignment1;
            alignment2 = '-' + alignment2;
            i--;
        } else if (traceback[i][j] === '←') {
            alignment1 = '-' + alignment1;
            alignment2 = seq2[j - 1] + alignment2;
            j--;
        } else {
            break;
        }
    }

    return `Alignment 1: ${alignment1}\nAlignment 2: ${alignment2}\nScore: ${dp[maxPos[0]][maxPos[1]]}`;
}

// Codon Usage Bias
function analyzeCodonUsage() {
    const sequence = document.getElementById('codonSequence').value.trim().toUpperCase();
    const codons = {};
    const standardCodons = { "AUG": "Methionine", "UUU": "Phenylalanine", "UUC": "Phenylalanine" }; // Example

    if (sequence) {
        for (let i = 0; i < sequence.length - 2; i += 3) {
            const codon = sequence.slice(i, i + 3);
            codons[codon] = (codons[codon] || 0) + 1;
        }

        let result = "Codon Usage:\n";
        for (const [codon, count] of Object.entries(codons)) {
            result += `${codon}: ${count} (Standard: ${standardCodons[codon] || "Unknown"})\n`;
        }

        document.getElementById('codonResult').innerText = result;
    } else {
        document.getElementById('codonResult').innerText = "Please enter a valid sequence.";
    }
}

// Motif Search
function searchMotif() {
    const sequence = document.getElementById('motifSequence').value.trim().toUpperCase();
    const motif = document.getElementById('motif').value.trim().toUpperCase();

    if (sequence && motif) {
        const indices = [];
        for (let i = 0; i <= sequence.length - motif.length; i++) {
            if (sequence.slice(i, i + motif.length) === motif) {
                indices.push(i + 1);
            }
        }

        document.getElementById('motifResult').innerHTML = indices.length
            ? `Motif found at positions: ${indices.join(', ')}`
            : "Motif not found.";
    } else {
        document.getElementById('motifResult').innerText = "Please enter a valid sequence and motif.";
    }
}
