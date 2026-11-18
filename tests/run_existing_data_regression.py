"""Optional regression against the user's original paired 95 Hz acceleration records."""
from pathlib import Path
import argparse
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import wave_packet as wp

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=wp.ROOT/'result/existing_data_regression')
    args=parser.parse_args()
    c=wp.load_config(wp.ROOT/'examples/existing_95Hz_acceleration.json')
    post=wp.import_file('packet_post_regression',wp.PROJECT/'postprocess_rt/postprocess_rt.py')
    target=args.data_root/'Planewave_upper_60_P_95Hz_addition';reference=args.data_root/'Planewave_upper_60_P_95Hz_uniform'
    result=post.analyze(c,target,reference,args.output,target/'sema',reference/'sema')
    expected={'R_P':(.5457738171027556,.5418582005646481),
              'R_SV':(.7222666931069892,.7192374723576512),
              'T_SV':(.354831991637492,.3556181413559274)}
    for row in result:
        np.testing.assert_allclose([row['packet_coefficient'],row['weighted_theory']],expected[row['phase']],rtol=2e-9,atol=1e-11)
    print('Existing-data regression passed: all three measured and weighted-theory coefficients reproduced.')
