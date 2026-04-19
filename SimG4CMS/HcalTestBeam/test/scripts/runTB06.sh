#!/bin/bash

python SimG4CMS/HcalTestBeam/test/python/run_tb06_all_cfg.py FTFP_BERT_EMM pi- $1 RR >& pi-.out &
python SimG4CMS/HcalTestBeam/test/python/run_tb06_kaon_cfg.py FTFP_BERT_EMM kaon+ $1 RR >& kaon+.out &
python SimG4CMS/HcalTestBeam/test/python/run_tb06_all_cfg.py FTFP_BERT_EMM pbar $1 RR >& pbar.out &
python SimG4CMS/HcalTestBeam/test/python/run_tb06_all_cfg.py FTFP_BERT_EMM p $1 RR >& p.out &
python SimG4CMS/HcalTestBeam/test/python/run_tb06_all_cfg.py FTFP_BERT_EMM pi+ $1 RR >& pi+.out &
python SimG4CMS/HcalTestBeam/test/python/run_tb06_kaon_cfg.py FTFP_BERT_EMM kaon- $1 RR >& kaon-.out &


#
