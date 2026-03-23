---
name: commit
description: Automatically generate a commit message from repository changes and create a git commit. Use when the user says /commit, asks to commit current changes, or requests automatic commit message generation.
---

# Commit

## Purpose
Create a safe, high-quality commit with an auto-generated message based on current repository changes.

## When to Use
- User explicitly asks to commit changes.
- User uses `/commit`.
- User asks to auto-generate a commit message and submit commit.

## Workflow
1. Collect context with:
   - `git status --short`
   - `git diff` (staged + unstaged)
   - `git log -n 10 --oneline`
2. Analyze change type and intent:
   - feature / fix / refactor / docs / test / chore
   - summarize **why** and major impact
3. Generate commit message:
   - first line: concise summary
   - optional body: key rationale and scope
4. Stage changes:
   - include relevant tracked and untracked files
   - exclude obvious secrets and local temp artifacts
5. Commit using HEREDOC message format.
6. Run `git status` and report result.

## Commit Message Style
- Follow existing project history style.
- Keep summary specific and actionable.
- Prefer intent-oriented wording over file-list wording.

Example:

```text
align summary merge order with java hashset semantics

Ensure deterministic merge traversal to eliminate read-assignment drift
against Java CIRI3 baseline.
```

## Guardrails
- Do not commit if there are no effective changes.
- Do not commit secret files (e.g., `.env`, credentials, private keys).
- Do not use destructive git operations.
- If commit hook fails, fix issues and create a new commit (do not silently skip checks).
