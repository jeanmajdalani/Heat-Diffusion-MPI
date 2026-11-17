# Project-specific LLM logging policy

For every AI-assisted interaction in this repository:
1. Record the date/time, the team member, AI tool/model, complete user prompt and complete assistant response.
2. Record the source files created/modified and line ranges where practical.
3. Append to `llm_log/LOG.md` without removing previous entries.
4. Commit AI-log updates with code changes.
5. For ChatGPT web conversations, manually export/copy the complete conversation into `llm_log/`; this instruction file cannot automatically capture web-chat messages.

The course requires all AI interactions, including rejected suggestions and prior planning/setup chats. The student must understand and independently verify all submitted work.
