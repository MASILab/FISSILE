# FISSILE

I've compiled `nnMpEvenRotationMatrix_v2` for better performance, but I've only generated the MEX file for macOS.

**Non-Mac users**: Run these commands in MATLAB to compile

```matlab
euler_type = coder.typeof(0, [Inf, 3]);
l_type = coder.typeof(0);
codegen nnMpEvenRotationMatrix_v2 -args {euler_type, l_type} -o nnMpEvenRotationMatrix_v2
```

