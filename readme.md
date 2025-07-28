# overview
Here is the structure of the code:
 - there is a class representing the compressor
 - The compressor is made of stages, so there is a class representing the stage
 - the stage is made of rotor and stator, and there is a class representing the velocity triangle of a single row of blades.

TLDR: There will be a class inside a class inside a class.

The reason for this is readability.

# IF YOU HAVE ANY QUESTIONS TEXT ME PLEASE


# we are pretty much making a J85-GE-5
 - 20,000 RPM
 - 8 stages + 1 IGV

 - https://media.defense.gov/2020/Apr/23/2002287288/-1/-1/0/200417-F-VV067-1152.JPG
 - https://media.defense.gov/2020/Apr/23/2002287289/-1/-1/0/200417-F-VV067-1153.JPG

# design
 - first stage gets free vortex:
   -  $C_w \cdot r = constant$ for C_w = 0 (First stage specifically)
 - Subsequent stages have options:
   - free vortex blading ($\lambda m$ = 0.5)
   - constant reaction blading ($\lambda m$ = 0.5) (What?)
   - Exponential blading ($\lambda m$ = 0.5) (What???)
   - to do: reference notes
 - investegate for all three different options
 - Free vortex design - on third stage rotor: 
   - get root mean and tip radii
   - get root mean and tip velocities
   - get root mean and tip betas - (beta is the incidence angle of the air on the rotor blade with reference to axial - arctangent(blade v / axial v))
   - free vortex blading condition - $C_{w2r} \cdot r = constant$
   - remember when you calculated the $C_{w2m} (no)
   - annulus area shrinks, so the calculate the angle the passage walls make with axial.