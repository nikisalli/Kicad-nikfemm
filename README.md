This kicad plugin uses nikfemm simulator to plot power density and voltage on PCB geometries.
Please note that this is just a demo and I tested it on linux systems only (but you could manage to make it work on windows quite easily).
Also you'll need to install all the dependencies by yourself.

# installation
to install the nikfemm python module:
- ```git clone https://github.com/nikisalli/nikfemm.git```
- ```cd nikfemm/python```
- ```pip install .```
(you'll need cmake, make and a compiler installed on your system to build the wheel)

to install the dependencies for the plugin
- ```git clone https://github.com/nikisalli/Kicad-nikfemm```
- ```cd Kicad-nikfemm```
- ```pip install -r requirements.txt```

to install the plugin:
- go to the kicad plugin manager
- click the install from file button
- select the zip in the main folder of this repo



https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/08204902-4b45-483e-9cc3-d83dded151fe



https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/887e64f4-ccd9-4500-953b-6109ba01c733

![photo_2024-02-17_17-42-46](https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/5f928cd1-b6d8-43fb-bb2b-1d909fd0ca02)
![photo_2024-02-17_17-42-50](https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/09856a6f-e10e-440d-a62d-2a2e25d1758f)
![photo_2024-02-17_17-42-53](https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/40370596-4531-4a54-bde8-1c819008af45)
![photo_2024-02-17_17-42-55](https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/a8f5245f-cf2b-4013-8580-aaecd93df261)
![photo_2024-02-17_17-42-42](https://github.com/nikisalli/Kicad-nikfemm/assets/31286021/519f64f9-bdfd-4e6d-bf7e-2b164c8c2c3f)
