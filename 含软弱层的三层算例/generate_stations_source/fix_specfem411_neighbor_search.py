#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python 3.10：预览/修复SPECFEM3D邻单元搜索索引。默认只预览，--apply备份后应用。

仅修改 src/specfem3D/locate_point.f90 的 find_best_neighbor 子程序。
不修改 mesh、DATABASES_MPI、FORCESOLUTION、STATIONS 或 Par_file。
"""
from pathlib import Path
from datetime import datetime
import argparse
import difflib
import os
import re
import stat
import tempfile

# ================= 用户编辑区 =================
SPECFEM_ROOT = '/share/SPECFEM3D_LYH/specfem3d-4.1.1'
# 默认只预览；实际应用时使用命令行 --apply。应用前自动保存原文件。
# ================= 用户编辑区结束 =================


def patched_text(text):
    newline='\r\n' if '\r\n' in text else '\n'
    normalized=text.replace('\r\n','\n')
    pattern=r'(?ims)^  subroutine find_best_neighbor\(.*?^  end subroutine find_best_neighbor[ \t]*$'
    matches=list(re.finditer(pattern,normalized))
    if len(matches)!=1:
        raise ValueError('Cannot uniquely identify find_best_neighbor; no file changed.')
    match=matches[0];body=match.group(0)
    marker='neighbor_start = neighbors_xadj(ispec_selected)'
    if marker in body and 'ientry = neighbor_start + ii' in body:
        return text
    replacements=[
        ('  integer :: ii,ientry','  integer :: ii,ientry,neighbor_start'),
        ('  num_neighbors = neighbors_xadj(ispec_selected+1) - neighbors_xadj(ispec_selected)',
         '  ! Freeze this neighbor list: ispec_selected may change during the search.\n'
         '  neighbor_start = neighbors_xadj(ispec_selected)\n'
         '  num_neighbors = neighbors_xadj(ispec_selected+1) - neighbor_start'),
        ('    ientry = neighbors_xadj(ispec_selected) + ii','    ientry = neighbor_start + ii'),
    ]
    for old,new in replacements:
        if body.count(old+'\n')!=1:
            raise ValueError('Source differs from the verified 4.1.1 pattern; no file changed: '+old.strip())
        body=body.replace(old+'\n',new+'\n',1)
    updated=normalized[:match.start()]+body+normalized[match.end():]
    return updated.replace('\n',newline)


def apply_to_file(path,apply=False):
    path=Path(path).resolve(strict=True)
    raw=path.read_bytes();old=raw.decode('utf-8');new=patched_text(old)
    if old==new:
        print('Already fixed:',path);return False
    print(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),
                                    fromfile=str(path),tofile=str(path)+' (fixed)')),end='')
    if not apply:
        print('\nPreview only. Use --apply to back up and apply this change.');return False
    if path.read_bytes()!=raw:
        raise RuntimeError('Source changed during inspection; nothing written.')
    backup=path.with_name(path.name+'.before_neighbor_fix_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    with backup.open('xb') as f:f.write(raw)
    mode=stat.S_IMODE(path.stat().st_mode)
    temp_path=None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix='.locate_point_fix_',delete=False) as f:
            temp_path=Path(f.name);f.write(new.encode('utf-8'))
        os.chmod(temp_path,mode)
        os.replace(temp_path,path)
    finally:
        if temp_path is not None and temp_path.exists():temp_path.unlink()
    print('\nApplied:',path,'\nBackup:',backup)
    print('Rebuild xspecfem3D before running. This script does not rebuild or run the solver.')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(SPECFEM_ROOT))
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    try:apply_to_file(args.root/'src/specfem3D/locate_point.f90',args.apply)
    except (OSError,ValueError,RuntimeError) as exc:parser.exit(1,str(exc)+'\n')
