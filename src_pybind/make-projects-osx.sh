cd core
mkdir buildXcode 
cd buildXcode
cmake -G Xcode -DPYTHON_EXECUTABLE:PATH=../../../../myenv/bin/python3  ..

cd gl
mkdir buildXcode 
cd buildXcode
cmake -G Xcode -DPYTHON_EXECUTABLE:PATH=../../../../myenv/bin/python3  ..

cd eigen
mkdir buildXcode 
cd buildXcode
cmake -G Xcode -DPYTHON_EXECUTABLE:PATH=../../../../myenv/bin/python3  ..