#!/bin/sh
graphlan_annotate.py MiceGraphlan.txt MiceGraphlan.xml --annot MiceGraphlan_anno.txt
graphlan.py MiceGraphlan.xml MiceGraphlan.pdf --dpi 300 --size 17 --pad 0

graphlan_annotate.py Stoma_CC_IBD_Control.txt Stoma_CC_IBD_Control.xml --annot Stoma_CC_IBD_Control_anno.txt
graphlan.py Stoma_CC_IBD_Control.xml Stoma_CC_IBD_Control.pdf --dpi 300 --size 8 --pad 0
