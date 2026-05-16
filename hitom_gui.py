"""
Hi-TOM GUI Launcher

A simple tkinter wrapper around hitom_analyze.py.  It lets you pick input
files, set the abundance threshold, and run the analysis with live log
output — no command-line required.

Usage:
    python hitom_gui.py
"""

import os
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

# Import the analysis engine
import hitom_analyze


class HiTomGUI:
    def __init__(self, root):
        self.root = root
        root.title("Hi-TOM Mutation Detection")
        root.geometry("720x560")
        root.resizable(False, False)

        # ---- Styling ----
        style = ttk.Style(root)
        style.configure("TButton", padding=4, font=("Segoe UI", 10))
        style.configure("TLabel", font=("Segoe UI", 10))
        style.configure("TEntry", font=("Segoe UI", 10))

        # ---- Input Files ----
        frame = ttk.LabelFrame(root, text="Input Files", padding=10)
        frame.pack(fill="x", padx=10, pady=5)

        # R1
        ttk.Label(frame, text="R1 FASTQ (.fq.gz):").grid(row=0, column=0, sticky="w")
        self.r1_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.r1_var, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse…", command=lambda: self._pick_fq(self.r1_var)).grid(row=0, column=2)

        # R2
        ttk.Label(frame, text="R2 FASTQ (.fq.gz):").grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.r2_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.r2_var, width=60).grid(row=1, column=1, padx=5, pady=(8, 0))
        ttk.Button(frame, text="Browse…", command=lambda: self._pick_fq(self.r2_var)).grid(row=1, column=2, pady=(8, 0))

        # Reference
        ttk.Label(frame, text="Reference FASTA:").grid(row=2, column=0, sticky="w", pady=(8, 0))
        self.ref_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.ref_var, width=60).grid(row=2, column=1, padx=5, pady=(8, 0))
        ttk.Button(frame, text="Browse…", command=self._pick_ref).grid(row=2, column=2, pady=(8, 0))

        # BAM
        ttk.Label(frame, text="Sorted BAM:").grid(row=3, column=0, sticky="w", pady=(8, 0))
        self.bam_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.bam_var, width=60).grid(row=3, column=1, padx=5, pady=(8, 0))
        ttk.Button(frame, text="Browse…", command=self._pick_bam).grid(row=3, column=2, pady=(8, 0))

        # ---- Parameters ----
        param_frame = ttk.LabelFrame(root, text="Parameters", padding=10)
        param_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(param_frame, text="Abundance threshold (%) :").grid(row=0, column=0, sticky="w")
        self.thresh_var = tk.StringVar(value="5.0")
        ttk.Entry(param_frame, textvariable=self.thresh_var, width=8).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Label(param_frame, text="Variants below this percentage are filtered out.").grid(row=0, column=2, sticky="w")

        # ---- Output ----
        out_frame = ttk.LabelFrame(root, text="Output", padding=10)
        out_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(out_frame, text="Output prefix:").grid(row=0, column=0, sticky="w")
        self.prefix_var = tk.StringVar(value="results")
        ttk.Entry(out_frame, textvariable=self.prefix_var, width=40).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Label(out_frame, text="→ <prefix>-sequence.tsv & <prefix>-genotype.tsv").grid(row=0, column=2, sticky="w")

        # ---- Run button ----
        btn_frame = ttk.Frame(root)
        btn_frame.pack(fill="x", padx=10, pady=5)
        self.run_btn = ttk.Button(btn_frame, text="▶ Run Analysis", command=self._run)
        self.run_btn.pack(side="left")
        self.progress = ttk.Progressbar(btn_frame, mode="indeterminate", length=200)
        self.progress.pack(side="right", padx=10)

        # ---- Log window ----
        log_frame = ttk.LabelFrame(root, text="Log", padding=5)
        log_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.log = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, state="disabled",
                                             font=("Consolas", 9))
        self.log.pack(fill="both", expand=True)

    # ------------------------------------------------------------------
    # File pickers
    # ------------------------------------------------------------------
    def _pick_fq(self, var):
        path = filedialog.askopenfilename(
            title="Select FASTQ file", filetypes=[("FASTQ", "*.fq.gz *.fastq.gz"), ("All files", "*.*")]
        )
        if path:
            var.set(path)
            # Auto-fill R2 if R1 is chosen and looks like paired-end
            if var is self.r1_var and "_1" in path:
                r2 = path.replace("_1", "_2")
                if os.path.exists(r2):
                    self.r2_var.set(r2)

    def _pick_ref(self):
        path = filedialog.askopenfilename(
            title="Select reference FASTA", filetypes=[("FASTA", "*.fasta *.fa"), ("All files", "*.*")]
        )
        if path:
            self.ref_var.set(path)

    def _pick_bam(self):
        path = filedialog.askopenfilename(
            title="Select sorted BAM", filetypes=[("BAM", "*.bam"), ("All files", "*.*")]
        )
        if path:
            self.bam_var.set(path)
            # Auto-fill reference if same directory has a .fasta
            d = os.path.dirname(path)
            for f in os.listdir(d):
                if f.endswith("_reference_sequence.fasta"):
                    self.ref_var.set(os.path.join(d, f))
                    break

    # ------------------------------------------------------------------
    # Logging helper
    # ------------------------------------------------------------------
    def _log(self, msg):
        self.log.configure(state="normal")
        self.log.insert(tk.END, msg + "\n")
        self.log.see(tk.END)
        self.log.configure(state="disabled")

    # ------------------------------------------------------------------
    # Run analysis in a background thread
    # ------------------------------------------------------------------
    def _run(self):
        r1 = self.r1_var.get().strip()
        r2 = self.r2_var.get().strip()
        ref = self.ref_var.get().strip()
        bam = self.bam_var.get().strip()
        prefix = self.prefix_var.get().strip()

        if not all([r1, r2, ref, bam]):
            messagebox.showerror("Missing input", "Please select all input files.")
            return

        if not os.path.isfile(r1):
            messagebox.showerror("File not found", f"R1 file not found:\n{r1}")
            return
        if not os.path.isfile(r2):
            messagebox.showerror("File not found", f"R2 file not found:\n{r2}")
            return
        if not os.path.isfile(ref):
            messagebox.showerror("File not found", f"Reference file not found:\n{ref}")
            return
        if not os.path.isfile(bam):
            messagebox.showerror("File not found", f"BAM file not found:\n{bam}")
            return
        bai = bam + ".bai"
        if not os.path.isfile(bai):
            messagebox.showwarning("Index missing", f"BAM index not found:\n{bai}\n\n"
                                                    f"Please run: samtools index {bam}")
            return

        try:
            thresh = float(self.thresh_var.get().strip()) / 100.0
            if not (0 < thresh < 1):
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid threshold", "Abundance threshold must be between 0.1 and 99.9.")
            return

        out_seq = f"{prefix}-sequence.tsv"
        out_geno = f"{prefix}-genotype.tsv"

        self.run_btn.configure(state="disabled")
        self.progress.start()
        self._log("=" * 50)
        self._log("Starting Hi-TOM analysis…")

        def worker():
            try:
                hitom_analyze.process_dataset(
                    dataset_num=1,
                    r1_path=r1,
                    r2_path=r2,
                    bam_path=bam,
                    ref_path=ref,
                    out_seq=out_seq,
                    out_geno=out_geno,
                    abundance_threshold=thresh,
                    log_callback=self._log,
                )
                self.root.after(0, self._on_success, out_seq, out_geno)
            except Exception as e:
                self.root.after(0, self._on_error, str(e))

        threading.Thread(target=worker, daemon=True).start()

    def _on_success(self, out_seq, out_geno):
        self.progress.stop()
        self.run_btn.configure(state="normal")
        messagebox.showinfo("Analysis complete",
                            f"Outputs written:\n\n{out_seq}\n{out_geno}")

    def _on_error(self, msg):
        self.progress.stop()
        self.run_btn.configure(state="normal")
        messagebox.showerror("Analysis failed", msg)


def main():
    root = tk.Tk()
    gui = HiTomGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
