## FBX Reader for GMS2
**This project was not working in IDE version v2023.11.1.129 (Steam), probably due to GMS new file format. I've convert fbx extensions into Gamemaker script instead, it should now works in newer GM version.**

This is example of skeletal animation for Game Maker Studio 2 readed from fbx using gml.
![screenshot](https://i.imgur.com/vqBnRhg.jpg)

Unfortunatelly not all animations can be loaded using this readed. Please check this list of points below. If all of the points fits your fbx file then more chanses it can be readed using this asset.
* fbx version 7400 or higher
* binary file
* one mesh in file (not strict)
* one animation in file (not strict)
* animation only skeletal (not strict)

I tested this asset only using fbx 7400 version exported from Blender.  
**Look out for any error when exporting fbx using Blender, sometime it fail to calculate tangent or something, this script will fail when reading undefined value.**

### Credits
* Original script by **Alexander Kondyrev** ([twitter.com/alsekond](https://twitter.com/alsekond))
* All fbx examples was downloaded from: [sketchfab.com](https://sketchfab.com)
