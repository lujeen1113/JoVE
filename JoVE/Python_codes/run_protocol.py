import sys
from BeautifiedNitin2 import setSysVar,SpeedSlideWindow
#setSysVar()
x=input('total sensor node amount is :')
total_Amount=int(x)
sectIncrAmount=(total_Amount/6)-6
sectIncrAmount=int(sectIncrAmount)
print(sectIncrAmount)
SpeedSlideWindow(2,0.008,sectIncrAmount)

