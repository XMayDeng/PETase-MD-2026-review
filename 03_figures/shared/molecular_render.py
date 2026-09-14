"""Execute only declared local molecular scenes, without cached-image alternatives."""
import os
from pathlib import Path
import re
import subprocess

EXPECTED_VERSION='1.13.dev202606262139'


def render_scenes(root, record, executable, work):
    binary=Path(executable).resolve()
    if not binary.is_file():raise FileNotFoundError('Explicit ChimeraX executable does not exist')
    env=os.environ.copy()
    env.update(XDG_CONFIG_HOME=str(work/'config'),XDG_CACHE_HOME=str(work/'cache'),
               XDG_DATA_HOME=str(work/'data'),PETASE_RENDER_WORK=str(work),PYTHONDONTWRITEBYTECODE='1')
    version=subprocess.run([str(binary),'--version'],capture_output=True,text=True,check=True,env=env,timeout=60)
    if f'UCSF ChimeraX version: {EXPECTED_VERSION} ' not in version.stdout+version.stderr:
        raise RuntimeError('ChimeraX version differs from the tested scene-rendering environment')
    for index,job in enumerate(record['render_scenes']):
        source=root/job['scene']
        commands=source.read_text()
        def resolve(match):
            base=root if match.group(1)=='PACKAGE' else work
            path=base/match.group(2)
            if '"' in str(path):raise ValueError('Double quotes in paths are not supported by the scene format')
            return '"'+str(path)+'"'
        commands=re.sub(r'@(PACKAGE|WORK)@/([^\s"]+)',resolve,commands)
        if '@PACKAGE@' in commands or '@WORK@' in commands:raise ValueError('Unresolved scene path token')
        script=work/f'scene_{index}.cxc';script.write_text(commands)
        result=subprocess.run([str(binary),'--nogui','--offscreen','--usedefaults','--notools',str(script)],
            cwd=work,env=env,capture_output=True,text=True,timeout=900)
        log=result.stdout+result.stderr
        (work/f'scene_{index}.log').write_text(log)
        if result.returncode or '\nERROR:' in log or 'Traceback (most recent call last)' in log:
            raise RuntimeError(f'Molecular scene failed ({job["scene"]}):\n{log[-6000:]}')
        for name in job['outputs']:
            path=work/name
            if not path.is_file() or path.stat().st_size==0:
                raise RuntimeError('Required scene output missing: '+name)
