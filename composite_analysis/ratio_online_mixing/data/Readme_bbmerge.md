# README: Running `merge_all_fastq.py` with WSL and BBMerge

## **📌 Overview**

This script automates the process of merging paired-end FASTQ files using **BBMerge** from the **BBMap** suite inside WSL (Windows Subsystem for Linux). It:

1. **Finds all **``** FASTQ files** in the dataset folder.
2. **Identifies the corresponding **``** files**.
3. **Runs **`` to merge the paired-end reads.
4. **Saves merged and unassembled files** in the same directory.

---

## **1️⃣ Requirements**

Before running the script, ensure you have the following installed:

### **✅ Windows Subsystem for Linux (WSL)**

1. Install WSL (Ubuntu recommended):

   ```powershell
   wsl --install
   ```

2. Restart your computer.

3. Open **Windows Terminal** and enter WSL:

   ```bash
   wsl
   ```

4. **Update Ubuntu**:

   ```bash
   sudo apt update && sudo apt upgrade -y
   ```

---

### **✅ Install BBMap in WSL**

1. Download BBMap from the official source:\
   [https://sourceforge.net/projects/bbmap/](https://sourceforge.net/projects/bbmap/)
2. Extract it into **Windows** under:
   ```
   C:\Users\User\Downloads\BBMap_39.10
   ```
3. Now, access it from WSL:
   ```bash
   cd /mnt/c/Users/User/Downloads/BBMap_39.10/bbmap
   ls -l
   ```
   You should see `bbmerge.sh` listed.

---

### **✅ Ensure Your FASTQ Files Are in the Correct Location**

Your **FASTQ files should be in Windows**, inside:

```
C:\Users\User\PycharmProjects\dna_storage_comb_syn_analysis\composite_analysis\ratio_online_mixing\data
```

WSL accesses this as:

```
/mnt/c/Users/User/PycharmProjects/dna_storage_comb_syn_analysis/composite_analysis/ratio_online_mixing/data
```

---

## **2️⃣ How to Run the Script**

### **🔹 Step 1: Open Terminal in PyCharm**

1. Open **PyCharm**.
2. Navigate to the **Terminal** tab (bottom panel).

### **🔹 Step 2: Run the Script**

```bash
python merge_all_fastq.py
```

The script will: ✔ Find `_R1_` and `_R2_` files\
✔ Verify their existence\
✔ Run `bbmerge.sh` for each file pair\
✔ Save the merged output in the same directory

---

## **3️⃣ Expected Output**

```
✅ Found R1 files:
  - /mnt/c/.../24587_S1_L001_R1_001.fastq.gz
  - /mnt/c/.../24588_S2_L001_R1_001.fastq.gz
...
🔎 Checking for R2 file: /mnt/c/.../24587_S1_L001_R2_001.fastq.gz
   ⏳ Running: wsl ls /mnt/c/.../24587_S1_L001_R2_001.fastq.gz
   📜 Output: /mnt/c/.../24587_S1_L001_R2_001.fastq.gz
   🔢 Return Code: 0
✅ R2 file exists!
🚀 Running bbmerge.sh for 24587_S1_L001...
   🔧 Command: wsl bash /mnt/c/.../bbmerge.sh in1=... in2=... out=...
```

After all files are processed, the merged files appear in:

```
C:\Users\User\PycharmProjects\dna_storage_comb_syn_analysis\composite_analysis\ratio_online_mixing\data
```

---

## **4️⃣ Debugging Issues**

If the script **skips files** or **fails**, follow these steps:

### **🛠️ Check If WSL Can See Your Files**

Run:

```bash
wsl ls -l /mnt/c/Users/User/PycharmProjects/dna_storage_comb_syn_analysis/composite_analysis/ratio_online_mixing/data
```

If `_R1_` or `_R2_` files **do not appear**, they may be in the wrong location.

---

### **🛠️ Manually Run **``** for Debugging**

If `bbmerge.sh` does not run, test it manually:

```bash
wsl bash /mnt/c/Users/User/Downloads/BBMap_39.10/bbmap/bbmerge.sh \
  in1=/mnt/c/Users/User/PycharmProjects/dna_storage_comb_syn_analysis/composite_analysis/ratio_online_mixing/data/24587_S1_L001_R1_001.fastq.gz \
  in2=/mnt/c/Users/User/PycharmProjects/dna_storage_comb_syn_analysis/composite_analysis/ratio_online_mixing/data/24587_S1_L001_R2_001.fastq.gz \
  out=/mnt/c/Users/User/PycharmProjects/dna_storage_comb_syn_analysis/composite_analysis/ratio_online_mixing/data/24587_S1_L001_merged.fastq.gz
```

If this **fails**, check if BBMap is installed correctly.

---

### **🛠️ Restart WSL (If Needed)**

```powershell
wsl --shutdown
wsl
```

Then **try running the script again**.

---

## **5️⃣ Summary of Key Commands**

| **Task**                  | **Command**                                       |
| ------------------------- | ------------------------------------------------- |
| Open WSL                  | `wsl`                                             |
| Update WSL Ubuntu         | `sudo apt update && sudo apt upgrade -y`          |
| List FASTQ files          | `ls -l /mnt/c/.../data/`                          |
| Test if R2 file exists    | `wsl ls /mnt/c/.../24587_S1_L001_R2_001.fastq.gz` |
| Run `bbmerge.sh` manually | `wsl bash /mnt/c/.../bbmerge.sh ...`              |
| Restart WSL               | `wsl --shutdown`                                  |

---

## **🎯 Final Notes**

✔ **Keep files in Windows** (`C:\Users\User\...`).\
✔ **Access them through WSL** (`/mnt/c/Users/User/...`).\
✔ **Use **``** to check files** before running the script.\
✔ **Run **``** manually** if needed.

This README ensures you can **re-run the script anytime without issues**. 🚀

