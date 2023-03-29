pro readsdssmet

fits_read, 'gal_fiboh_dr7_v5_2.fits', data 

;IDL> print, data
; 194 199 204 205 194 199 204 205 194 199 204 205 194 199 204 205 194 199 204 205
; 194 199 204 205 194 199 204 205 192  88 249 153 160   0   0   0
 

; IDL> fits_read, 'gal_fiboh_dr7_v5_2.fits', data, /noscale
; IDL> help, data
; DATA            BYTE      = Array[36, 927552]


;; contains BYTE array corresponding to Median, P16, P84, P2.5, P97.5, mode, avg, and entropy

; IDL> print, data[0:3,0:5]
print, uint(data[0:3,0:5])  
print, float(data[0:3,0:5],4,6)
print, float(data[0:3,0:5],0,6)
f=-99.9
b=byte(f,0,4,1)
print, b
print, float(data[3:0,0:5],0,6)
reverse=[3,2,1,0]
print, float(data[reverse,0:5],0,6)
print, float(data[reverse,*],0,6)  
print, float(data[reverse,*],0,100)
median=float(data[reverse,*],0,92700)
reverse=[7,6,5,4]                    
p16=float(data[reverse,*],0,92700)   
reverse=[11,10,9,8]               
p84=float(data[reverse,*],0,92700)
print((p84[where(p84 ne -99.9)]-p16[where(p84 ne -99.9)])/2.0)

print, ((p84[where(p84 ne -99.9)]-p16[where(p84 ne -99.9)])/2.0)
print, median((p84[where(p84 ne -99.9)]-p16[where(p84 ne -99.9)])/2.0)

print, max((p84[where(p84 ne -99.9)]-p16[where(p84 ne -99.9)])/2.0)   

plothist, ((p84[where(p84 ne -99.9)]-p16[where(p84 ne -99.9)])/2.0), bin=0.01


stop
end