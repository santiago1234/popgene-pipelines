# Gnomix code

Get the most up to date version

```bash
git clone https://github.com/AI-sandbox/gnomix.git
```

## Notes on `bin`


- I changed this line from [gnomix](https://github.com/AI-sandbox/gnomix/blob/2eca7b0a9ddf67c569cf95088b764c61b1050932/src/postprocess.py#L186C30-L186C41) to this: ` "chm": chm` because if not it fails my chromosome names have the format `chrN2`
- Also, I changed [another line](https://github.com/AI-sandbox/gnomix/blob/2eca7b0a9ddf67c569cf95088b764c61b1050932/src/postprocess.py#L207C47-L207C79) to `sample_file_name = os.path.join(root, sample.replace(".", "_").strip() + ".bed")`

