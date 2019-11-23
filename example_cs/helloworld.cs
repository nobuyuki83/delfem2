using System;
using System.Runtime.InteropServices;

class Program
{
  private static class NativeMethods
  {
    [DllImport("../src_dll/build/libdelfem2.so", EntryPoint="add")]
    public static extern int add(int a, int b);
  }

  static void Main(string[] args)
  {
    int a = NativeMethods.add(2,3);
    Console.WriteLine(a);
  }
}	
