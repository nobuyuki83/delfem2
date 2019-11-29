using System;
using System.Runtime.InteropServices;

class Program
{

#if (__APPLE__)
　private const string DLL_PATH = "../src_dll/build/libdelfem2.dylib";  
#else
　private const string DLL_PATH = "../src_dll/build/libdelfem2.so";
#endif

  private static class NativeMethods
  {   

    [DllImport(DLL_PATH, EntryPoint="add")]
    public static extern int add(int a, int b);
  }

  static void Main(string[] args)
  {
    int a = NativeMethods.add(2,3);
    Console.WriteLine(a);
    Console.WriteLine(Environment.OSVersion.VersionString);
    Console.WriteLine(DLL_PATH);
  }
}	
