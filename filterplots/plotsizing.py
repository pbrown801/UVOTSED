#method for sizing plots correctly

#matplotlib takes plot size in inches, so use this function to convert
#input is in pts desired, or you can include an if statement with a certain name like journal to have that size built in
def setFigSizeInches(width, fraction=1, subplots=(1,1)):
    #can include presets like this so you don't have to type in exact pts every time
    if(width == 'journal'):
        width_pt = 242.26653
        heightfactor = 1.5 #makes the height look nicer imo, than pure gr
    else:
        width_pt = width
        heightfactor = 1    
    print("starting")
    figwidthpt = width_pt * fraction
    #to convert from pt to inches
    ippt = 1 / 72.27 
    #golden ratio for asthetics
    gr = (5**.5 -1)/2
    #width and height in inches 
    figwidthin = figwidthpt * ippt 
    figheightin = figwidthin * gr * heightfactor * (subplots[0] / subplots[1])
    #make sure to use plt.tight_layout() to keep labels nice
    figdim = (figwidthin, figheightin)
    return figdim

#use this function if for some reason you don't need inch conversion
def setFigSize(width, fraction=1):
    print("starting")
    figwidthpt = width * fraction
    #to convert from pt to inches
    #golden ratio for asthetics
    gr = (5**.5 -1)/2
    #width and height in inches 
    figwidthin =  242.26653 #figwidthpt
    figheightin = figwidthin * gr

    figdim = (figwidthin, figheightin)
    return figdim