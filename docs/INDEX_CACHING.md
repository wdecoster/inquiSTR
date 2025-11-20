# Index Caching for Remote Files

## Overview

When accessing remote BAM/CRAM files (via HTTP/HTTPS/FTP/S3), inquiSTR automatically manages index file caching to improve performance and reduce redundant downloads.

## How It Works

### 1. Local Index Check (Priority #1)

Before accessing a remote file, inquiSTR checks if an index file already exists **in the current directory**:

- For `file.bam`: checks for `file.bam.bai` or `file.bai`
- For `file.cram`: checks for `file.cram.crai` or `file.crai`

If a local index is found, it will be used instead of downloading from the remote location.

**Example:**

```bash
# If you have the index locally, just place it in the working directory
wget https://example.com/file.cram.crai
inquiSTR call https://example.com/file.cram --reference genome.fa ...
# inquiSTR will use the local file.cram.crai
```

### 2. Remote Index Caching (Automatic)

If no local index is found and you're accessing a remote file, inquiSTR automatically sets up a cache directory for htslib:

- **Default location**: `~/.cache/inquistr/` (Linux/Mac) or `.inquistr_cache/` (if HOME is not set)
- **Behavior**: Once downloaded, the index is cached and reused in subsequent runs
- **Manual override**: Set `HTS_CACHE_LOCATION` environment variable before running inquiSTR

**Example:**

```bash
# First run: downloads and caches the index
inquiSTR unmapped https://example.com/file.cram --reference genome.fa

# Second run: uses cached index (much faster!)
inquiSTR unmapped https://example.com/file.cram --reference genome.fa
```

### 3. Manual Cache Location

You can specify a custom cache directory by setting the `HTS_CACHE_LOCATION` environment variable:

```bash
export HTS_CACHE_LOCATION=/scratch/my_cache
inquiSTR call https://example.com/file.cram --reference genome.fa ...
```

## Performance Benefits

- **First access**: Index is downloaded once and cached
- **Subsequent access**: Index is read from cache (instant)
- **Multiple files**: Each index is cached separately, no re-downloads needed
- **Automatic cleanup**: Files older than 30 days are automatically removed on each run

## Cache Management

### Automatic Cleanup

inquiSTR automatically removes cached index files older than 30 days every time it accesses a remote file. This prevents the cache from growing indefinitely.

**Customizing cleanup age:**

```bash
# Keep cached files for 60 days instead of 30
export INQUISTR_CACHE_MAX_AGE_DAYS=60
inquiSTR call https://example.com/file.cram --reference genome.fa ...
```

### Manual Cache Management

**View cache contents:**

```bash
ls -lh ~/.cache/inquistr/
```

**Clean cache manually (dry run):**

```bash
inquiSTR clean-cache --dry-run
```

**Clean files older than 7 days:**

```bash
inquiSTR clean-cache --max-age-days 7
```

**Clean entire cache:**

```bash
inquiSTR clean-cache --all
# or manually:
rm -rf ~/.cache/inquistr/
```

### Disabling Cache

For one-time analyses where you don't want to cache the index:

```bash
export INQUISTR_NO_CACHE=1
inquiSTR call https://example.com/file.cram --reference genome.fa ...
```

This is useful when:

- Processing a file you'll never use again
- Testing different files frequently
- Working with temporary/ephemeral data

## Technical Details

- inquiSTR checks for local indexes using common extensions (`.bai`, `.crai`, `.bam.bai`, `.cram.crai`)
- Cache directory is created automatically if it doesn't exist
- The caching mechanism is provided by htslib (the underlying C library)
- Cache is shared across all htslib-based tools that use the same `HTS_CACHE_LOCATION`

## Troubleshooting

### Cache Not Working?

Check if `HTS_CACHE_LOCATION` is already set in your environment:

```bash
echo $HTS_CACHE_LOCATION
```

If it's set to a location you don't have write access to, unset it or change it:

```bash
unset HTS_CACHE_LOCATION
# or
export HTS_CACHE_LOCATION=/path/to/writable/directory
```

### Index Download Errors?

If the remote index cannot be downloaded, inquiSTR will fall back to:

1. Using a local index if available
2. Sequential reading (without index) if possible
3. Error message if index is required for the operation

## Best Practices

1. **For frequent access**: Let inquiSTR cache the index automatically (default behavior)
2. **For one-time use**: Either download the index manually to current directory, or set `INQUISTR_NO_CACHE=1`
3. **For shared systems**: Set a shared cache location with `HTS_CACHE_LOCATION`
4. **For cache management**: Use `inquiSTR clean-cache --dry-run` to see what will be cleaned
5. **For debugging**: Use `RUST_LOG=debug` to see cache operations:

```bash
RUST_LOG=debug inquiSTR unmapped https://example.com/file.cram --reference genome.fa 2>&1 | grep -i cache
```

## Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `INQUISTR_NO_CACHE` | Disable caching entirely | Not set (caching enabled) |
| `INQUISTR_CACHE_MAX_AGE_DAYS` | Days to keep cached files | 30 |
| `HTS_CACHE_LOCATION` | Override cache directory | `~/.cache/inquistr/` |

## Examples

**Normal usage (automatic caching):**

```bash
inquiSTR call https://example.com/file.cram -R repeats.bed --reference genome.fa
```

**Disable caching for one-time analysis:**

```bash
INQUISTR_NO_CACHE=1 inquiSTR call https://example.com/file.cram -R repeats.bed --reference genome.fa
```

**Keep cache files for 90 days:**

```bash
export INQUISTR_CACHE_MAX_AGE_DAYS=90
inquiSTR call https://example.com/file.cram -R repeats.bed --reference genome.fa
```

**Check cache status:**

```bash
inquiSTR clean-cache --dry-run
```
