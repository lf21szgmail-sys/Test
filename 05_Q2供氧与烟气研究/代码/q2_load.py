
"""只读提取哈希已核验的单表 XLSB。遇到公式、错误或未知单元格类型时拒绝读取。"""
import csv
import hashlib
import struct
import zipfile
from pathlib import Path

COLUMNS = (0, 2048, 2049, 2050)
HERE = Path(__file__).resolve().parent

def records(stream):
    while first := stream.read(1):
        kind = first[0] & 127
        if first[0] & 128:
            kind |= stream.read(1)[0] << 7
        size = 0
        for shift in range(0, 28, 7):
            digit = stream.read(1)[0]
            size |= (digit & 127) << shift
            if digit < 128:
                break
        yield kind, stream.read(size)

def read_book(path, heat):
    with zipfile.ZipFile(path) as book:
        from xml.etree import ElementTree as ET
        from posixpath import normpath, join
        rels = ET.fromstring(book.read("xl/_rels/workbook.bin.rels"))
        sheets = [r.attrib["Target"] for r in rels if r.attrib["Type"].endswith("/worksheet")]
        if len(sheets) != 1:
            raise ValueError(f"{heat}: expected one worksheet")
        target = sheets[0]
        sheet = target.lstrip("/") if target.startswith("/") else normpath(join("xl", target))
        names = []
        with book.open("xl/workbook.bin") as stream:
            for kind, data in records(stream):
                if kind == 156:  # BrtBundleSh, after state and sheet ID
                    length = struct.unpack_from("<I", data, 8)[0]
                    offset = 12 + length * 2
                    count = struct.unpack_from("<I", data, offset)[0]
                    names.append(data[offset+4:offset+4+count*2].decode("utf-16le"))
        if names != [heat]:
            raise ValueError((path.name, names, heat))
        strings = []
        if "xl/sharedStrings.bin" in book.namelist():
            with book.open("xl/sharedStrings.bin") as stream:
                for kind, data in records(stream):
                    if kind == 19:
                        count = struct.unpack_from("<I", data, 1)[0]
                        strings.append(data[5:5+count*2].decode("utf-16le"))
        rows, header = [], {}
        row_index = -1
        with book.open(sheet) as stream:
            for kind, data in records(stream):
                if kind == 0:
                    index = struct.unpack_from("<I", data)[0]
                    if index != row_index + 1:
                        raise ValueError(f"{heat}: nonconsecutive Excel rows")
                    row_index = index
                    if index:
                        rows.append([None] * 4)
                elif 1 <= kind <= 11:
                    col = struct.unpack_from("<I", data)[0]
                    if col not in COLUMNS and row_index != 0:
                        continue
                    if kind == 1:
                        value = None
                    elif kind == 2:
                        rk = struct.unpack_from("<i", data, 8)[0]
                        value = rk >> 2 if rk & 2 else struct.unpack("<d", struct.pack("<Ii", 0, rk & ~3))[0]
                        if rk & 1:
                            value /= 100
                    elif kind == 5:
                        value = struct.unpack_from("<d", data, 8)[0]
                    elif kind == 6:
                        count = struct.unpack_from("<I", data, 8)[0]
                        value = data[12:12+count*2].decode("utf-16le")
                    elif kind == 7:
                        value = strings[struct.unpack_from("<I", data, 8)[0]]
                    else:
                        raise ValueError((heat, row_index, col, kind))
                    if row_index == 0:
                        if col in header:
                            raise ValueError((heat, col, "duplicate header"))
                        header[col] = value
                    else:
                        rows[-1][COLUMNS.index(col)] = value
        if set(header) != set(range(2051)):
            raise ValueError((heat, len(header)))
        labels = [header[c] for c in COLUMNS]
        cleaned = [str(s).lower().replace("₂","2").replace(" ","") for s in labels]
        if cleaned[0] not in {"sequence","序列","序列sequence"} or cleaned[1] not in {"累积耗氧量","yqxh","cumulativeoxygenconsumption"} or cleaned[2] not in {"co2浓度","co2concentration"} or cleaned[3] not in {"co浓度","coconcentration"}:
            raise ValueError((heat, labels))
        return labels, rows

def manifest(data_dir):
    import re
    with (HERE / "版本清单.csv").open(encoding="utf-8-sig", newline="") as stream:
        records_ = list(csv.DictReader(stream))
    seen = set()
    for i, item in enumerate(records_, 1):
        name = item["原始文件名"]
        if int(item["附件编号"]) != i or (i == 1 and name != "附件1.炉号-6B03205.xlsb"):
            raise ValueError("Invalid manifest")
        if i > 1 and not name.startswith(f"Attachment {i}. "):
            raise ValueError(name)
        heat = re.search(r"(\d[A-Z]\d{5})\.xlsb$", name)[1]
        path = data_dir / name
        with path.open("rb") as stream:
            sha = hashlib.file_digest(stream, "sha256").hexdigest()
        if path.stat().st_size != int(item["文件字节数"]) or sha != item["SHA256"] or heat in seen:
            raise ValueError((path.name, sha))
        seen.add(heat)
        yield i, heat, path, sha, int(item["采样点数"])
    if len(seen) != 10:
        raise ValueError("Expected ten distinct heats")

if __name__ == "__main__":
    import sys
    for i, heat, path, sha, n in list(manifest(Path(sys.argv[1]))):
        if i in (1,2,5):
            labels, rows = read_book(path, heat)
            print(i, heat, len(rows), labels, rows[0], rows[-1])
